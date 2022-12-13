import logging
import pandas as pd
import numpy as np
from typing import Union
import arcpy
from arcpy.sa import Raster, Filter, Con, CreateConstantRaster, MajorityFilter
from file_functions import table_to_csv, delete_gis_files, cmd, \
    err_info, spatial_license
from create_centerline import make_centerline
from create_station_lines import create_station_lines_function
import os
import shutil


@err_info
@spatial_license
def lidar_footprint(
    lasbin: str,
    lidardir: str,
    spatialref_shp: str,
) -> str:
    """This function converts LAZ files to LAS file format as well as producing a LiDAR extent polygon.
    in_folder must be a directory containing nothing but raw LAZ files
    spatial_ref must be an ArcGIS spatial reference object with units of feet.
    las_tools_bin must be the location of the 'bin' folder installed with LAStools by rapidlasso
    Returns: A shapefile w/ LiDAR coverage to be used to make a ground polygon for LAStools processing"""
    files_in_direct = [f for f in os.listdir(
        lidardir) if os.path.isfile(os.path.join(lidardir, f))]

    laspath = lidardir + '\\las_files'

    if not os.path.exists(laspath):
        os.makedirs(laspath)

    # Initiate temp files folder formatted for LAStools
    temp_files = lidardir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    in_spatial_ref = arcpy.Describe(spatialref_shp).spatialReference

    # Convert laz files to LAS files
    for f in files_in_direct:
        if f[-4:] == ".laz":

            # Correct format, can alter between browse() input and default
            if lasbin[-1] != 'n':
                lasbin = lasbin[:-1]

            cmd("%s\\laszip.exe -i %s\\%s -o %s\\%s_noprj.las" %
                (lasbin, lidardir, f, laspath, f[:-4]))
            logging.info("%s\\laszip.exe -i %s\\%s -o %s\\%s_noprj.las" %
                         (lasbin, lidardir, f, laspath, f[:-4]))
            cmd("%s\\las2las.exe -i %s\\%s_noprj.las -o %s\\%s.las" %
                (lasbin, laspath, f[:-4], laspath, f[:-4]))
            logging.info("%s\\las2las.exe -i %s\\%s_noprj.las -o %s\\%s.las" %
                         (lasbin, laspath, f[:-4], laspath, f[:-4]))

    files_in_laspath = [f for f in os.listdir(
        laspath) if os.path.isfile(os.path.join(laspath, f))]

    # Delete unnecessary index files
    for f in files_in_laspath:
        if f[-4:] == 'lasx':
            os.remove(laspath + "\\%s" % f)

        if f[-5] == 'j':
            os.remove(laspath + "\\%s" % f)

    raw_las_dataset = arcpy.CreateLasDataset_management(
        laspath,
        lidardir + "\\raw_las_dataset.lasd",
        spatial_reference=in_spatial_ref,
        compute_stats=True,
    )
    lidar_ras = CreateConstantRaster(1, extent=raw_las_dataset)
    lidar_footprint = arcpy.RasterToPolygon_conversion(
        lidar_ras,
        lidardir + '\\las_footprint.shp',
    )

    return lidar_footprint


@err_info
@spatial_license
def define_ground_polygon(
    lidar_footprint: str,
    lidardir: str,
    naipdir: str,
    ndvi_thresh: float,
    aoi_shp: str,
) -> str:
    """This function takes the defined lidar footprint from the lidar_footprint() function, as well as a defined NAIP imagery location (in .jpg2)
    and makes a polygon of vegeation using a NDVI threshold of >0.4. This polygon is erased from the lidar footprint to give a ground_polygon used
    to define processing settings"""

    # Set processing extent to the LiDAR data extent
    arcpy.env.extent = lidar_footprint
    in_spatial_ref = arcpy.Describe(lidar_footprint).spatialReference

    # Find NAIP imagery in folder
    naip_imagery = [f for f in os.listdir(
        naipdir) if os.path.isfile(os.path.join(naipdir, f))]

    # Initiate temp files folder
    temp_files = lidardir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    if len(naip_imagery) > 1:
        add_to_mosaic = [naipdir + "\\" + f for f in naip_imagery]
        naip_imagery = arcpy.MosaicToNewRaster_management(
            add_to_mosaic,
            output_location=lidardir,
            raster_dataset_name_with_extension="NAIP_mos.tif",
            coordinate_system_for_the_raster=in_spatial_ref,
            number_of_bands=4,
        )
    else:
        naip_imagery = (naipdir + "\\%s" % naip_imagery[0])
        naip_imagery = arcpy.ProjectRaster_management(
            naip_imagery,
            lidardir + "\\NAIP_prj.tif",
            in_spatial_ref,
        )

    # Extract bands 1 (red) and 4 (NIR)
    red_lyr = arcpy.MakeRasterLayer_management(
        naip_imagery,
        temp_files + "\\rd_lyr",
        band_index=0,
    )
    nir_lyr = arcpy.MakeRasterLayer_management(
        naip_imagery,
        temp_files + "\\nr_lyr",
        band_index=4,
    )

    red_lyr = arcpy.SaveToLayerFile_management(
        red_lyr,
        temp_files + "\\red_ras.lyr",
    )
    nir_lyr = arcpy.SaveToLayerFile_management(
        nir_lyr,
        temp_files + "\\nir_ras.lyr",
    )

    red_ras = arcpy.CopyRaster_management(
        red_lyr,
        temp_files + "\\red_ras.tif",
        format="TIFF",
    )
    nir_ras = arcpy.CopyRaster_management(
        nir_lyr,
        temp_files + "\\nir_ras.tif",
        format="TIFF",
    )

    red_ras = Raster(red_ras)
    nir_ras = Raster(nir_ras)

    # Calculate ndvi and generate polygon delineating values > ndvi_thresh
    ndvi = lidardir + "\\NDVI.tif"
    ndvi_ras = ((nir_ras - red_ras) / (nir_ras + red_ras))
    ndvi_ras.save(ndvi)

    veg_ras_raw = Con(Raster(ndvi) >= ndvi_thresh, 1)
    veg_ras_raw.save(temp_files + "\\veg_ras_raw.tif")

    veg_ras = MajorityFilter(
        veg_ras_raw,
        "EIGHT",
        "MAJORITY",
    )
    veg_ras.save(temp_files + "\\veg_ras.tif")

    veg_poly = arcpy.RasterToPolygon_conversion(
        veg_ras,
        lidardir + "\\veg_poly_ndvi.shp",
        simplify="FALSE",
    )

    # Make polygon representing bare ground
    if aoi_shp != '':
        ground_poly = arcpy.Erase_analysis(
            lidar_footprint,
            veg_poly,
            temp_files + "\\ground_poly_full.shp",
        )
        aoi_prj = arcpy.Project_management(
            aoi_shp,
            temp_files + "\\aoi_prj_to_inref.shp",
            out_coor_system=in_spatial_ref,
        )
        ground_poly = arcpy.Clip_analysis(
            ground_poly,
            aoi_prj,
            lidardir + "\\ground_poly.shp",
        )

    else:
        ground_poly = arcpy.Erase_analysis(
            lidar_footprint,
            veg_poly,
            lidardir + "\\ground_poly.shp",
        )

    ground_poly = arcpy.DefineProjection_management(
        ground_poly,
        in_spatial_ref,
    )
    logging.info("AOI bare-ground polygon @ %s" % ground_poly)

    return ground_poly


@err_info
@spatial_license
def lidar_to_raster(
    lidardir: str,
    spatialref_shp: str,
    aoi_shp: str,
    sample_meth: str,
    tri_meth: str,
    void_meth: str,
    m_cell_size: Union[float, int] = 1,
) -> str:
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
        logging.info('LAS units are Meters')

    elif in_spatial_ref.linearUnitName == 'Foot_US':
        cell_size = (3.28 * m_cell_size)
        logging.info('LAS units are Feet')

    else:
        return logging.info('Linear unit name for %s uncertain, please use a PROJECTED COORDINATE SYSTEM' % os.path.basename(in_spatial_ref))

    # Set up interpolation method string
    if sample_meth == 'BINNING':
        method_str = '%s AVERAGE %s' % (sample_meth, void_meth)

    else:
        method_str = "%s %s NO_THINNING MAXIMUM 0" % (sample_meth, tri_meth)
    logging.info('Methods: %s' % method_str)

    no_prj_dem = temp_files + '\\noprj_dem.tif'
    las_dataset = arcpy.CreateLasDataset_management(
        ground_lasdir,
        out_las,
        spatial_reference=in_spatial_ref,
        compute_stats=True,
    )
    lidar_raster = arcpy.LasDatasetToRaster_conversion(
        las_dataset,
        value_field='ELEVATION',
        data_type='FLOAT',
        interpolation_type=method_str,
        sampling_type='CELLSIZE',
        sampling_value=cell_size,
    )
    arcpy.CopyRaster_management(lidar_raster, no_prj_dem)
    arcpy.ProjectRaster_management(
        no_prj_dem,
        out_raster=out_dem,
        out_coor_system=out_spatial_ref,
    )

    logging.info("LAS -> DEM output @ %s" % out_dem)

    # Notify the user which units the DEM are in
    if out_spatial_ref.linearUnitName == 'Meter':
        logging.info('DEM units are Meters')

    elif out_spatial_ref.linearUnitName == 'Foot_US':
        logging.info('DEM units are Feet')

    else:
        logging.info('Linear unit name for %s uncertain, please use a PROJECTED COORDINATE SYSTEM' %
                     os.path.basename(out_spatial_ref))

    return out_dem


@err_info
@spatial_license
def detrend_prep(
    dem: str,
    flow_poly: str,
    aoi_shp: str,
    filt_passes: int,
    smooth_dist: Union[int, float],
    m_spacing: Union[int, float] = 1,
    centerline_verified: bool = False,
) -> str:
    """This function takes the Lidar raster, creates a least-cost thalweg centerline from a smoothed raster. Station points are
    generated along the centerline at defined spacing (1/20th of channel width is a starting point) which are given the values of the lidar raster.

    Args: raster_name, upstream flow polygon, spatial extent (can be raster), station point spacing in ft (3ft is default).
    Run first with centerline_verified=False and visually inspect. Run again w/ True to return the [station_points, elevation_table]"""

    # Set up environment and output folder
    spatial_ref = arcpy.Describe(aoi_shp).spatialReference
    arcpy.env.extent = dem
    dem_dir = os.path.dirname(dem)

    # Initiate temp files folder
    temp_files = dem_dir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    # Define input parameters
    params = [m_spacing,
              smooth_dist]  # First item defines XS length and spacing, second item described smoothing distance

    if not spatial_ref.linearUnitName == 'Meter':
        params = [int(i * 3) for i in params]

    filt_passes = int(filt_passes)

    if not centerline_verified:
        logging.info('Generating smooth thalweg centerline...')
        logging.info("Smoothing DEM w/ %sx low pass filters..." % filt_passes)
        ticker = 0
        filter_out = Filter(dem, "LOW")
        filter_out.save(temp_files + "\\filter_out%s" % ticker)

        while ticker < filt_passes:  # Apply an iterative low pass filter 15x to the raster to smooth the topography
            filter_out = Filter(
                (temp_files + "\\filter_out%s" % ticker),
                "LOW",
            )
            filter_out.save(temp_files + "\\filter_out%s" % (ticker + 1))
            ticker += 1
        smooth_ras = (dem_dir + "\\filt_ras.tif")
        filter_out.save(dem_dir + "\\filt_ras.tif")

        # Create least cost centerline from 15x filtered raster
        logging.info(
            "Smoothed DEM made, least-cost centerline being calculated...")
        lidar_foot = dem_dir + '\\las_footprint.shp'

        # check for LiDAR Footprint file
        if not os.path.exists(lidar_foot):
            logging.error(f'Could not for the previously generated las_footprint.shp file in {dem_dir} \
            ...please re-make the shapefile or move it back to the default folder.')
            logging.error(
                'WARNING: The following process may run but will propduce incorrect outputs!')

        make_centerline(
            smooth_ras,
            aoi_shp,
            lidar_foot,
            flow_poly,
            smooth_distance=10,
        )

        # Delete intermediate filtered rasters
        for ticker in range(filt_passes + 1):
            file = (temp_files + "\\filter_out%s" % ticker)
            if os.path.exists(file):
                try:
                    shutil.rmtree(file)
                except Exception:
                    logging.info("Could not remove %s " % file)
            else:
                logging.info("Path %s does not exist and can't be deleted...")
        logging.info('Done')

    else:
        logging.info('Generating thalweg elevation profile...')
        centerline = dem_dir + "\\thalweg_centerline.shp"

        # Define location of intermediate files, some of which will be deleted
        intermediates = [
            'thalweg_centerline_XS.shp',
            'thalweg_station_points.shp',
            'thalweg_station_points1.shp',
            'sp_elevation_table.dbf',
        ]
        intermediates = [temp_files + '\\%s' % i for i in intermediates]

        # Create a station point shapefile evenly sampling the thalweg centerline
        station_lines = create_station_lines_function(
            centerline,
            spacing=params[0],
            xs_length=params[0],
        )
        station_points = arcpy.Intersect_analysis(
            [intermediates[0], centerline],
            out_feature_class=intermediates[2],
            join_attributes="ALL",
            output_type="POINT",
        )
        station_points = arcpy.MultipartToSinglepart_management(
            station_points,
            intermediates[1],
        )
        station_points = arcpy.AddXY_management(station_points)

        # Extract elevation values from each station point, and export to a .csv file
        elevation_table = arcpy.ExtractValuesToTable_ga(
            station_points,
            in_rasters=dem,
            out_table=intermediates[3],
        )
        station_points = arcpy.JoinField_management(
            station_points,
            in_field="ORIG_FID",
            join_table=elevation_table,
            join_field="SrcID_Feat",
            fields=["Value"],
        )

        # Add fields to override, but first adjust detrending functions
        elevation_table = dem_dir + '\\xyz_elevation_table.csv'
        elevation_table = table_to_csv(
            input_table=station_points,
            csv_filepath=elevation_table,
            fld_to_remove_override=[
                'FID_thal_1',
                'Id_1',
                'InLine_FID',
                'ORIG_FID'
            ],
            keep_fields=[],
        )
        elevation_df = pd.read_csv(elevation_table)

        # Flip rows if upside down
        max_loc = elevation_df['LOCATION'].max()
        elevation_df.sort_values('LOCATION', inplace=True)

        if elevation_df.iloc[0]['Value'] < elevation_df.iloc[-1]['Value']:
            loc_list = elevation_df.loc[:, ['LOCATION']].squeeze().to_list()
            loc_np = np.array([int(max_loc - i) for i in loc_list])
            elevation_df['LOCATION'] = loc_np
            elevation_df.sort_values('LOCATION', inplace=True)
        elevation_df.to_csv(elevation_table)

        # Delete extra files
        for j in intermediates[2:]:
            delete_gis_files(j)

        logging.info("Thalweg elevation profile (.csv) @ %s " %
                     str(elevation_table))
        logging.info('Done')

        return elevation_table
