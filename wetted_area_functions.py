import os
import numpy as np
from typing import Union, List
from matplotlib import pyplot as plt
import arcpy
from arcpy.sa import Raster, Con, BoundaryClean, MajorityFilter
from arcpy.da import SearchCursor
from file_functions import get_label_units, string_to_list
from create_centerline import remove_spurs


def float_keyz_format(
    z: float,
    n: int = 1,
) -> str:
    """This function takes a float key z argument and retrusn its equivalent formatted string.
    ex: 5.3 -> 5p3, or 10.0 -> 10p0
    If the n parameter (default=1) is altered, more digits past the decimal as converted"""
    z_str = ''
    if z >= 10.0 and isinstance(z, float):
        z_str = (str(z)[0:2] + 'p' + str(z)[3:3+n])
    elif z < 10.0 and isinstance(z, float):
        z_str = (str(z)[0] + 'p' + str(z)[2:3+n])
    elif isinstance(z, int):
        z_str = str(z) + 'p0'

    try:
        return z_str
    except z_str == '':
        print('Key z list parameters not valid. Please fill list with int or float.')


def prep_small_inc(
    detrended_dem: str,
    max_stage: Union[float, int],
) -> str:
    """IN: Folder containing detrended DEM ras_detren.tif, an stage interval length, a maximum flood stafe height.
    RETURNS: Directory with wetted polygons"""
    # Set up an out directory for wetted area polygons
    in_dir = os.path.dirname(detrended_dem)
    out_dir = in_dir + '\\wetted_polygons'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Set up a new folder storing intermediate rasters used in centerline generation and for viewing
    int_files = out_dir + '\\wetted_area_rasters'
    if not os.path.exists(int_files):
        os.makedirs(int_files)

    # Use detrended dem spatial reference to determine stage intervals
    u, unit, spatial_ref = get_label_units(detrended_dem)

    if unit == 'Meter':
        interval = 0.03
        n = 2
    else:
        interval = 0.1
        n = 1
    print('Units are %s' % unit)

    # Set up range
    stages = np.arange(0, float(max_stage + interval), interval)
    stages = np.around(stages, n)

    # Make wetted area polygons at 0.1ft / 0.03m intervals
    print('Making wetted polygons...')
    in_ras = Raster(detrended_dem)

    for inc in stages:
        inc_str = float_keyz_format(inc, n)
        temp_names = [int_files + '\\noval_%s%s.tif' %
                      (inc_str, u), int_files + '\\dt_clp_%s%s.tif' % (inc_str, u)]
        out_name = out_dir + '\\wetted_poly_%s%s.shp' % (inc_str, u)

        # Only generate polygons that have not already been generated
        if not os.path.exists(out_name):
            # Create intermediate rasters and detrended dems clipped at each wetted interval
            wetted_ras = Con(in_ras <= inc, 1)
            clip_ras = Con(in_ras <= inc, in_ras)
            wetted_ras.save(temp_names[0])
            clip_ras.save(temp_names[1])

            # Turn the wetted area raster into a polygon, delete intermediate rasters
            arcpy.RasterToPolygon_conversion(
                in_raster=wetted_ras,
                out_polygon_features=out_name,
                simplify=False,
            )

    return out_dir


def pdf_cdf_plotting(
    in_dir: str,
    out_folder: str,
    max_stage: Union[float, int],
) -> List[str]:
    """Doc string goes here
    Returns: A list containing the locations of the three generated wetted area plots"""
    print('Wetted area vs stage height analysis initiated...')

    # Make new folder to hold plots
    out_folder = out_folder + '\\flow_stage_plots'

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Find all wetted area polygons in their out folder
    wetted_areas = []
    wetted_polys = [in_dir + '\\%s' % f for f in os.listdir(
        in_dir) if f[:11] == 'wetted_poly' and f[-4:] == '.shp']

    # Set units based on the end of the wetted polygons name
    if wetted_polys[0][-5] == 'm':
        interval = 0.03
        u = 'm'
    else:
        interval = 0.1
        u = 'ft'

    # Calculate the wetted area of each wetted area polygon in their folder
    print('Calculating wetted areas...')
    for poly in wetted_polys:
        poly_area = 0
        for row in SearchCursor(poly, ["SHAPE@AREA"]):
            poly_area += float(row[0])
        wetted_areas.append(poly_area)

    # Clean up any errors and sort wetted areas from smallest to largest
    wetted_areas = [i for i in wetted_areas if i is not None]
    wetted_areas.sort()

    # Clip list down to only include wetted areas below the selected max stage
    stages = np.arange(0, max_stage + interval, interval)
    wetted_areas = wetted_areas[:len(stages)]

    # Calculate the change in wetted area between stages
    print('Calculating d(wetted area)...')
    d_area = []
    for count, area in enumerate(wetted_areas):
        if count == 0:
            d_area.append(area)
        else:
            d_area.append(float(area-wetted_areas[count-1]))

    # Plot stage height (x axis) vs wetted area (y axis)
    print('Plotting...')
    x1 = stages
    y1 = np.array(wetted_areas)
    title1 = (out_folder + '\\cumulative_area.png')
    plt.figure()
    plt.plot(x1, y1)
    plt.xlabel('Flood stage height (%s)' % u, fontsize='small')
    plt.ylabel('Wetted area (sq %s)' % u, fontsize='small')
    plt.title('Cumulative wetted area chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, max_stage + 1, step=1), fontsize='x-small')
    plt.yticks(fontsize='x-small')
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.savefig(title1, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    # Plot the derivative of the previous plot: PDF plot, shows d(wetted area)
    x2 = np.arange(interval, max_stage + interval, interval)
    y2 = np.array(d_area[1:])
    title2 = (out_folder + '\\pdf_plot.png')
    plt.figure()
    plt.plot(x2, y2)
    plt.xlabel('Flood stage height (%s)' % u, fontsize='small')
    plt.ylabel('Change in area (sq %s)' % u, fontsize='small')
    plt.title('PDF: d(wetted area) plot')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage+1), step=1), fontsize='x-small')
    plt.yticks(fontsize='x-small')
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.savefig(title2, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    # Plot the z vs wetted area, but with wetted area as the x axis to represent an mean cross-sectional geometry
    x3 = np.array(wetted_areas)
    y3 = stages
    title3 = out_folder + '\\mean_XS_plot.png'
    plt.figure()
    plt.plot(x3, y3)
    plt.xlabel('Wetted area (sq %s)' % u, fontsize='small')
    plt.ylabel('Flood stage height (%s)' % u, fontsize='small')
    plt.title('Represents mean cross-sectional geometry')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#666666', linestyle='-')
    plt.xlim(0, max(x3))
    plt.ylim(0, max_stage)
    plt.xticks(np.arange(0, int(max(x3)), step=round(
        max(x3) / 10)), fontsize='x-small')
    plt.yticks(np.arange(0, int(max(y3)), step=1), fontsize='x-small')
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.savefig(title3, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    return [title1, title2, title3]


def stage_centerlines(dem, zs, drafting=True):
    """Inputs: A folder containing key stage wetted area polygons (including intermediate file folder). Zs, a list
    containing N number of stage heights (floats) or a string with key xs separated by commas (ex: '0.2,0.7,2.6')"""
    # convert from string to list if necessary

    if type(zs) == str:
        zs = string_to_list(zs, format='float')

    # set up directories
    dem_dir = os.path.dirname(dem)

    if len(dem_dir) == 0:
        print('Error: Please select valid detrended DEM file')
        return

    out_dir = dem_dir + '\\centerlines'
    wetted_dir = dem_dir + '\\wetted_polygons\\wetted_area_rasters'
    temp_files = dem_dir + '\\temp_files'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    drafts = []

    # set up messages
    if drafting:
        messages = [
            'Generating draft center-lines for flow stage heights %s...' % zs,
            'Draft center-lines @ %s. \nManually edit (must see documentation) and run the next step. \nDone' % out_dir,
        ]
    else:
        messages = [
            'Generating final center-lines for flow stage heights %s...' % zs,
            'Final center-lines @ %s \nDone.' % out_dir,
        ]

    print(messages[0])

    # set up units string
    spatial_ref = arcpy.Describe(dem).spatialReference
    unit = spatial_ref.linearUnitName
    if unit == 'Meter':
        u = 'm'
        spur_lim = 100
        smooth = 20
    else:
        u = 'ft'
        spur_lim = (100 * 3.28)
        smooth = 60

    # majority filter, boundary clean, raster to polygon, polygon to centerline, remove spurs
    if drafting:
        for i, z in enumerate(zs):

            z_str = float_keyz_format(z)
            in_name = wetted_dir + '\\noval_%s%s.tif' % (z_str, u)
            out_name = out_dir + '\\%s%s_centerline_draft.shp' % (z_str, u)

            mf = MajorityFilter(in_name, 'EIGHT')
            bc = BoundaryClean(mf)

            temp_poly = temp_files + '\\sp%s.shp' % i  # smoothed polygon
            arcpy.RasterToPolygon_conversion(bc, temp_poly)

            w_spurs = temp_files + '\\%s%s_spur_cl.shp' % (z_str, u)
            rm_spur = w_spurs.replace('.shp', '_rm_spurs.shp')

            spurs = str(
                arcpy.PolygonToCenterline_topographic(
                    temp_poly,
                    w_spurs,
                ),
            )
            remove_spurs(spurs, spur_length=spur_lim)

            arcpy.CopyFeatures_management(rm_spur, out_name)
            drafts.append(out_name)

        print(
            'Please see centerline_info.txt in %s for information about editing centerlines')
        # print message, make .txt file in the centerline folder with instructions on editing

    elif not drafting:
        for z in zs:
            z_str = float_keyz_format(z)
            draft = out_dir + '\\%s%s_centerline_draft.shp' % (z_str, u)
            out_name = draft.replace('_draft.shp', '.shp')
            diss = temp_files + \
                os.path.basename(draft).replace('_draft.shp', 'diss.shp')

            # make into multipart, then slightly smooth
            arcpy.Dissolve_management(
                draft,
                diss,
                dissolve_field='ObjectID',
            )
            arcpy.SmoothLine_cartography(
                diss,
                out_name,
                'PAEK',
                smooth,
            )
            arcpy.AddField_management(
                out_name,
                'Id',
                'Short',
            )

    print(messages[1])
    return out_dir
