import arcpy
import numpy as np
from arcpy import env
from arcpy.sa import *
import os
import matplotlib
from matplotlib import pyplot as plt
import file_functions
from file_functions import *


def float_keyz_format(z):
    '''This function takes a float key z argument and retrusn its equivalent formatted string.
    ex: 5.3 -> 5p3, or 10.0 -> 10p0'''

    z_str = ''
    if z >= 10.0 and isinstance(z, float):
        z_str = (str(z)[0:2] + 'p' + str(z)[3:])
    elif z < 10.0 and isinstance(z, float):
        z_str = (str(z)[0] + 'p' + str(z)[2:])
    elif isinstance(z, int):
        z_str = str(z) + 'p0'

    try:
        return z_str
    except z_str == '':
        print('Key z list parameters not valid. Please fill list with int or float.')


def prep_small_inc(detrended_dem, max_stage):
    """IN: Folder containing detrended DEM ras_detren.tif, an stage interval length, a maximum flood stafe height.
    RETURNS: None. This function creates a folder containing wetted polygons for a 0.1ft increments as well as a clippped detrended DEM and contours."""
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
    spatial_ref = arcpy.Describe(detrended_dem).spatialReference
    unit = spatial_ref.linearUnitName
    if not unit == 'Meter':
        interval = 0.03
        u = 'm'
    else:
        interval = 0.1
        u = 'ft'
    print('TESTING. Units are %s' % unit)

    # Set up range
    stages = np.arange(0, max_stage + interval, float(interval))

    # Make wetted area polygons at 0.1ft / 0.03m intervals
    print('Making wetted polygons...')
    in_ras = arcpy.sa.Raster(detrended_dem)

    for inc in stages:
        inc_str = float_keyz_format(inc)
        temp_names = [int_files + '\\noval_%s%s.tif' % (inc_str, u), int_files + '\\dt_clp_%s%s.tif' % (inc_str, u)]
        out_name = out_dir + '\\wetted_poly_%s%s.shp' % (inc_str, u)

        # Create intermediate rasters and detrended dems clipped at each wetted interval
        wetted_ras = arcpy.sa.Con(in_ras <= inc, 1)
        clip_ras = arcpy.sa.Con(in_ras <= inc, in_ras)
        wetted_ras.save(temp_names[0])
        clip_ras.save(temp_names[1])

        # Turn the wetted area raster into a polygon, delete intermediate rasters
        arcpy.RasterToPolygon_conversion(in_raster=wetted_ras, out_polygon_features=out_name, simplify=False)
    print('Done')

    return out_dir


def pdf_cdf_plotting(in_dir, out_folder, max_stage):
    """Doc string goes here
    Returns: A list containing the locations of the three generated wetted area plots"""
    print('Wetted area vs stage height analysis initiated...')

    # Find all wetted area polygons in their out folder
    wetted_areas = []
    wetted_polys = [in_dir + '\\%s' % f for f in os.listdir(in_dir) if f[:11] == 'wetted_poly']

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
        for row in arcpy.da.SearchCursor(poly, ["SHAPE@AREA"]):
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
    title1 = (out_folder + '\\wetted_areas_plot_small_inc.png')
    plt.figure()
    plt.plot(x1, y1)
    plt.xlabel('Flood stage height (%s)' % u)
    plt.ylabel('Wetted area (sq %s)' % u)
    plt.title('Cumulative wetted area chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, max_stage + 1, step=1))
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
    plt.xlabel('Flood stage height (%s)' % u)
    plt.ylabel('Change in area (sq %s)' % u)
    plt.title('PDF: d(wetted area) plot')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage+1), step=1))
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
    plt.xlabel('Wetted area (sq %s)' % u)
    plt.ylabel('Flood stage height (%s)' % u)
    plt.title('Represents mean cross-sectional geometry')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#666666', linestyle='-')
    plt.xlim(0, max(x3))
    plt.ylim(0, max_stage)
    plt.xticks(np.arange(0, int(max(x3)), step=round(max(x3) / 10)))
    plt.yticks(np.arange(0, int(max(y3)), step=1))
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.savefig(title3, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')
    print('Done')

    return [title1, title2, title3]


