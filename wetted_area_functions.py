import arcpy
from arcpy import env
from arcpy.sa import *
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
