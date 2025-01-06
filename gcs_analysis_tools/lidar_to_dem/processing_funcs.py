import sys
import logging
import os
from typing import Union
from pathlib import Path

import arcpy
from arcpy.sa import Raster, Con, CreateConstantRaster, MajorityFilter

sys.path.append(str(Path(__file__).parent.parent))
from utils import (
    cmd,
    err_info, 
    spatial_license,
)


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
    lidar_footprint = lidardir + '\\las_footprint.shp'
    arcpy.RasterToPolygon_conversion(
        lidar_ras,
        lidar_footprint,
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


