"""Tests the functionality of the first GUI window: LiDAR to DEM."""
import sys
import pytest
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.parent))
from gcs_analysis_tools.lidar_to_dem.gui_funcs import (
    lidar_prep,
    dem_generation,
)

def test_lidar_prep(ROOT_DIR: Path, TEST_DATA_DIR: Path):
    """Tests the preparation of LiDAR data for conversion to DEM."""
    # test w/ all inputs
    lasbin = str(ROOT_DIR / "LAStools" / "bin")
    lidardir = str(TEST_DATA_DIR)
    spatial_shp = str(TEST_DATA_DIR / "lidar_prj_shp.shp")
    lidar_prep(
        lasbin=lasbin,
        lidardir=lidardir,
        spatial_shp=spatial_shp,
        naip_folder=lidardir,
        ndvi_thresh=0.4,
        aoi_shp="",
    )
    assert True
    assert False

    # test w/o aoi_shp
    lidar_prep(
        lasbin=lasbin,
        lidardir=lidardir,
        spatial_shp=spatial_shp,
        naip_folder="",
        ndvi_thresh=0.4,
    )

    assert True

def test_dem_generation():
    """Tests the LASTools conversion of LiDAR to DEM."""
    pytest.skip("Not implemented yet.")
