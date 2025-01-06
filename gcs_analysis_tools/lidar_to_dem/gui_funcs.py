import logging
from arcpy import HillShade_3d

from .processing_funcs import lidar_footprint, define_ground_polygon, lidar_to_raster
from .las_tools_funcs import get_bare_ground_las

def lidar_prep(
    lasbin: str,
    lidardir: str,
    spatial_shp: str,
    naip_folder: str,
    ndvi_thresh: float = 0.4,
    aoi_shp: str = '',
) -> None:
    """This function is ran by the LiDAR Data prep tab.
     Outputs: A 1m (or other resolution) DEM modeling bare ground LiDAR returns"""

    logging.info('Unzipping LAZ files...')
    foot: str = lidar_footprint(
        lasbin,
        lidardir,
        spatial_shp,
    )
    logging.info('Done')

    logging.info('Generating inputs for LiDAR processing...')
    define_ground_polygon(
        foot,
        lidardir,
        naip_folder,
        ndvi_thresh,
        aoi_shp,
    )
    logging.info('Done')


def dem_generation(
        lastoolsdir,
        lidardir,
        ground_poly,
        cores,
        units_code,
        keep_orig_pts,
        coarse_step,
        coarse_bulge,
        coarse_spike,
        coarse_down_spike,
        coarse_offset,
        fine_step,
        fine_bulge,
        fine_spike,
        fine_down_spike,
        fine_offset,
        aoi_shp,
        dem_resolution,
        dem_method,
        tri_meth,
        void_meth,
):
    """This function used LAStools to generate LAS file tiles prepresenitng bare ground, and then converts them
    to a high resolution DEM (1m resolution is default)"""

    # We carry input spatial ref over from the above process, but we should still convert from shp to ref object
    logging.info('Processing LiDAR to remove vegetation points...')
    las_folder = lidardir + '\\las_files\\'
    get_bare_ground_las(
        lastoolsdir + '\\',
        las_folder,
        ground_poly,
        cores,
        units_code,
        keep_orig_pts,
        coarse_step,
        coarse_bulge,
        coarse_spike,
        coarse_down_spike,
        coarse_offset,
        fine_step,
        fine_bulge,
        fine_spike,
        fine_down_spike,
        fine_offset,
    )
    logging.info('Done')

    logging.info('Generating a %sm resolution DEM...' % dem_resolution)
    dem = lidar_to_raster(
        lidardir,
        ground_poly,
        aoi_shp,
        dem_method,
        tri_meth,
        void_meth,
        m_cell_size=float(dem_resolution),
    )
    logging.info('Done')

    logging.info('Generating hillshade raster for the DEM...')
    hill_out = lidardir + '\\hillshade.tif'
    HillShade_3d(
        dem,
        hill_out,
    )
    logging.info('Done')


