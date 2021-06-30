import arcpy
from arcpy import env
from arcpy.sa import *
import file_functions
from file_functions import *
import create_centerline_GUI
import create_station_lines
from create_station_lines import create_station_lines_function
import os
from os import listdir
from os.path import isfile, join
import xlrd
import shutil
from openpyxl.workbook import Workbook
from openpyxl.reader.excel import load_workbook, InvalidFileException


def lidar_footptint(lasbin, lidardir, spatialref_shp):
    """This function converts LAZ files to LAS file format as well as producing a LiDAR extent polygon.
    in_folder must be a directory containing nothing but raw LAZ files
    spatial_ref must be an ArcGIS spatial reference object with units of feet.
    las_tools_bin must be the location of the 'bin' folder installed with LAStools by rapidlasso
    Returns: A shapefile w/ LiDAR coverage to be used to make a ground polygon for LAStools processing"""
    files_in_direct = [f for f in listdir(lidardir) if isfile(join(lidardir, f))]

    laspath = lidardir + '\\las_files'

    if not os.path.exists(laspath):
        os.makedirs(laspath)

    # Initiate temp files folder formatted for LAStools
    temp_files = lidardir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    spatial_ref = arcpy.Describe(spatialref_shp).spatialReference

    try:
        # Convert laz files to LAS files
        for f in files_in_direct:
            if f[-4:] == ".laz":

                # Correct format, can alter between browse() input and default
                if lasbin[-1] != 'n':
                    lasbin = lasbin[:-1]

                cmd("%s\\laszip.exe -i %s\\%s -o %s\\%s_noprj.las" % (lasbin, lidardir, f, laspath, f[:-4]))
                print("%s\\laszip.exe -i %s\\%s -o %s\\%s_noprj.las" % (lasbin, lidardir, f, laspath, f[:-4]))
                cmd("%s\\las2las.exe -i %s\\%s_noprj.las -o %s\\%s.las" % (lasbin, laspath, f[:-4], laspath, f[:-4]))
                print("%s\\las2las.exe -i %s\\%s_noprj.las -o %s\\%s.las" % (lasbin, laspath, f[:-4], laspath, f[:-4]))

        files_in_laspath = [f for f in listdir(laspath) if isfile(join(laspath, f))]

        # Delete unnecessary index files
        for f in files_in_laspath:
            if f[-4:] == 'lasx':
                os.remove(laspath + "\\%s" % f)

            if f[-5] == 'j':
                os.remove(laspath + "\\%s" % f)

        raw_las_dataset = arcpy.CreateLasDataset_management(laspath, lidardir + "\\raw_las_dataset.lasd",
                                                            spatial_reference=spatial_ref, compute_stats=True)
        lidar_footprint = arcpy.PointFileInformation_3d(raw_las_dataset, temp_files + "\\las_footprint_pre_dissolve",
                                                        "LAS", input_coordinate_system=spatial_ref)
        lidar_footprint = arcpy.Dissolve_management(lidar_footprint, lidardir + "\\las_footprint")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

    return lidar_footprint


def define_ground_polygon(lidar_footprint, lidardir, spatialref_shp, naipdir, ndvi_thresh, aoi_shp):
    """This function takes the defined lidar footprint from the lidar_footprint() function, as well as a defined NAIP imagery location (in .jpg2)
    and makes a polygon of vegeation using a NDVI threshold of >0.4. This polygon is erased from the lidar footprint to give a ground_polygon used
    to define processing settings"""

    # Set processing extent to the LiDAR data extent
    arcpy.env.extent = lidar_footprint
    spatial_ref = arcpy.Describe(aoi_shp).spatialReference

    # Find NAIP imagery in folder
    naip_imagery = [f for f in listdir(naipdir) if isfile(join(naipdir, f))]

    # Initiate temp files folder
    temp_files = lidardir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    if len(naip_imagery) > 1:
        add_to_mosaic = [naipdir + "\\" + f for f in naip_imagery]
        naip_imagery = arcpy.MosaicToNewRaster_management(add_to_mosaic, output_location=lidardir,
                                                          raster_dataset_name_with_extension="NAIP_mos.tif",
                                                          number_of_bands=4)
    else:
        naip_imagery = (naipdir + "\\%s" % naip_imagery[0])

    try:
        # Project the NAIP data and extract bands 1(red) and 4(NIR)
        naip_imagery = arcpy.ProjectRaster_management(naip_imagery, lidardir + "\\NAIP_prj.tif", spatial_ref)
        red_lyr = arcpy.MakeRasterLayer_management(naip_imagery, temp_files + "\\rd_lyr", band_index=0)
        nir_lyr = arcpy.MakeRasterLayer_management(naip_imagery, temp_files + "\\nr_lyr", band_index=4)

        red_lyr = arcpy.SaveToLayerFile_management(red_lyr, temp_files + "\\red_ras.lyr")
        nir_lyr = arcpy.SaveToLayerFile_management(nir_lyr, temp_files + "\\nir_ras.lyr")

        red_ras = arcpy.CopyRaster_management(red_lyr, temp_files + "\\red_ras.tif", format="TIFF")
        nir_ras = arcpy.CopyRaster_management(nir_lyr, temp_files + "\\nir_ras.tif", format="TIFF")

        red_ras = Raster(red_ras)
        nir_ras = Raster(nir_ras)

        # Calculate ndvi and generate polygon delineating values > ndvi_thresh
        ndvi = ((nir_ras - red_ras) / (nir_ras + red_ras))
        ndvi.save(lidardir + "//NDVI.tif")

        veg_ras_raw = Con(arcpy.sa.Raster(ndvi) >= ndvi_thresh, 1)
        veg_ras_raw.save(temp_files + "//veg_ras_raw.tif")
        veg_ras = MajorityFilter(veg_ras_raw, "EIGHT", "MAJORITY")
        veg_ras.save(temp_files + "//veg_ras.tif")
        veg_poly = arcpy.RasterToPolygon_conversion(veg_ras, lidardir + "//veg_poly_ndvi.shp", simplify=FALSE)

        # Make polygon representing bare ground
        if aoi_shp != '':
            ground_poly = arcpy.Erase_analysis(lidar_footprint, veg_poly, temp_files + "//ground_poly_full.shp")
            centerline_buff_prj = arcpy.Project_management(aoi_shp, temp_files + "//aoi_prj_to_inref.shp",
                                                           out_coor_system=spatial_ref)
            ground_poly = arcpy.Clip_analysis(ground_poly, aoi_shp, lidardir + "//ground_poly.shp")

        else:
            ground_poly = arcpy.Erase_analysis(lidar_footprint, veg_poly, lidardir + "//ground_poly.shp")

        print("AOI bare-ground polygon @ %s" % ground_poly)

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())


def lidar_to_raster(lidardir, spatialref_shp, aoi_shp, sample_meth, tri_meth, void_meth,  m_cell_size=1):
    """Converts processed LAS files to a LAS dataset, and then to a raster with cell size of 1m
    Args: Folder containing LAS files, desired cell size in meters (default is 1m), and ft spatial reference
    Returns: Raster name for use in detrending """
    # Create variables with relevant folders
    lasdir = lidardir + '\\las_files'
    ground_lasdir = lasdir + '\\09_ground_rm_duplicates'

    # Create addresses for generated .lasd, .tiff files
    out_dem = lidardir + "\\las_dem.tif"
    out_las = lasdir + '\\las_dataset.lasd'

    # Initiate temp files folder
    temp_files = lidardir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    # Set up output spatial reference and convert units if necessary
    in_spatial_ref = arcpy.Describe(spatialref_shp).spatialReference
    out_spatial_ref = arcpy.Describe(aoi_shp).spatialReference

    if in_spatial_ref.linearUnitName == 'Meter':
        cell_size = m_cell_size

    else:
        cell_size = (3.28 * m_cell_size)

    # Set up interpolation method string
    if sample_meth == 'BINNING':
        method_str = '%s AVERAGE %s' % (sample_meth, void_meth)

    else:
        method_str = "%s %s NO_THINNING MAXIMUM 0" % (sample_meth, tri_meth)
    print('Methods: %s' % method_str)

    try:
        no_prj_dem = temp_files + '\\noprj_dem.tif'
        las_dataset = arcpy.CreateLasDataset_management(ground_lasdir, out_las, spatial_reference=in_spatial_ref,
                                                        compute_stats=True)
        lidar_raster = arcpy.LasDatasetToRaster_conversion(las_dataset, value_field='ELEVATION', data_type='FLOAT',
                                                           interpolation_type=method_str, sampling_type='CELLSIZE', sampling_value=cell_size)
        tiff_lidar_raster = arcpy.CopyRaster_management(lidar_raster, out_dem)
        tiff_lidar_raster = arcpy.ProjectRaster_management(lidar_raster, out_raster=out_dem,
                                                           out_coor_system=out_spatial_ref)

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())
    print("LAS -> DEM output @ %s" % out_dem)

    return out_dem


def detrend_prep(raster_name, flow_polygon, spatial_extent, ft_spatial_ref, ft_spacing=3, use_filtered_ras=False,
                 centerline_verified=False):
    """This function takes the Lidar raster, creates a least-cost thalweg centerline from a smoothed raster. Station points are
    generated along the centerline at defined spacing (1/20th of channel width is a starting point) which are given the values of the lidar raster.

    Args: raster_name, upstream flow polygon, spatial extent (can be raster), station point spacing in ft (3ft is default).
    Run first with centerline_verified=False and visually inspect. Run again w/ True to return the [station_points, elevation_table]"""

    # arcpy.env.extent = spatial_extent
    raster_folder = os.path.dirname(raster_name)
    spacing = int(ft_spacing)
    xs_length = 5
    smooth_distance = 20
    filter_steps = 15

    if centerline_verified == False:
        print("Smoothing raster w/ 15x low pass filters...")
        ticker = 0
        filter_out = arcpy.sa.Filter(raster_name, "LOW")
        filter_out.save(raster_folder + "\\filter_out%s" % ticker)
        while ticker < filter_steps:  # Apply an iterative low pass filter 15x to the raster to smooth the topography
            filter_out = arcpy.sa.Filter((raster_folder + "\\filter_out%s" % ticker), "LOW")
            filter_out.save(raster_folder + "\\filter_out%s" % (ticker + 1))
            ticker += 1
        smooth_ras = (raster_folder + "\\filt_ras.tif")
        filter_out.save(raster_folder + "\\filt_ras.tif")

        print("Smoothed raster made, least-cost centerline being calculated...")
        least_cost_cl = create_centerline_GUI.least_cost_centerline(smooth_ras,
                                                                    upstream_source_poly)  # Create least cost centerline from 10x filtered raster
        least_cost_cl = create_centerline_GUI.remove_spurs(least_cost_cl, spur_length=10)
        centerline = create_centerline_GUI.smooth_centerline(least_cost_cl, smooth_distance=smooth_distance)

        for ticker in range(filter_steps + 1):  # Delete intermediate filtered rasters
            file = (raster_folder + "\\filter_out%s" % ticker)
            if os.path.exists(file):
                try:
                    shutil.rmtree(file)
                except:
                    print("Could not remove %s " % file)
            else:
                print("Path %s does not exist and can't be deleted...")
        print(
            "Intermediate files deleted, please manually verify centerline quality and define reach range... Call this function again w/ centerline_verified=True.")

    else:
        direct = os.path.dirname(flow_polygon)
        centerline = direct + "\\las_files\\centerline\\smooth_centerline.shp"
        station_lines = create_station_lines.create_station_lines_function(centerline, spacing=spacing,
                                                                           xs_length=xs_length, stage=[])

        station_lines = direct + ("\\las_files\\centerline\\smooth_centerline_XS_%sx%sft.shp" % (spacing, xs_length))
        print("Station lines file at: " + str(station_lines))
        print("Centerline at: " + str(centerline))

        station_points = arcpy.Intersect_analysis([station_lines, centerline], out_feature_class=(
                    direct + "\\station_points_%s_smooth_%s_spaced.shp" % (smooth_distance, spacing)),
                                                  join_attributes="ALL", output_type="POINT")
        station_points = arcpy.MultipartToSinglepart_management(station_points, (
                    direct + "\\raw_station_points_%s_smooth %s_spaced.shp" % (smooth_distance, spacing)))
        station_points = arcpy.AddXY_management(station_points)
        station_points = arcpy.Sort_management(station_points, out_dataset=(
                    direct + "\\XYZ_station_points_%s_smooth %s_spaced.shp" % (smooth_distance, spacing)),
                                               sort_field=[["LOCATION", "Ascending"]])

        if use_filtered_ras == True:
            raster_name = (raster_folder + "\\filt_ras.tif")
        elevation_table = arcpy.ExtractValuesToTable_ga(station_points, in_rasters=raster_name, out_table=(
                    direct + "\\sp_elevation_table_%s_smooth %s_spaced.dbf" % (smooth_distance, spacing)))
        station_points = arcpy.JoinField_management(station_points, in_field="ORIG_FID", join_table=elevation_table,
                                                    join_field="SrcID_Feat", fields=["Value"])

        elevation_table = direct + '\\XYZ_elevation_table.csv'
        elevation_table = file_functions.tableToCSV(input_table=station_points, csv_filepath=elevation_table,
                                                    fld_to_remove_override=[])  # Add fields to override, but first adjust detrending functions 'POINT_M', 'FID_smooth', 'Id', 'FID_smoo_1', 'InLine_FID', 'ORIG_FID'
        elevation_df = pd.read_csv(elevation_table)

        max_loc = elevation_df['LOCATION'].max()  # See if this thing works, hasn't been tested yet
        if elevation_df.iloc[0]['Value'] < elevation_df.iloc[-1]['Value']:
            loc_list = elevation_df.loc[:, [('LOCATION')]].squeeze().to_list()
            loc_np = np.array([int(max_loc - i) for i in loc_list])
            elevation_df['LOCATION'] = loc_np
            elevation_df.sort_values('LOCATION', inplace=True)
        elevation_df.to_csv(elevation_table)

        print("Station points shapefile at: " + str(station_points))
        print("Elevation table at: " + str(elevation_table))

        return [station_points, elevation_table]
