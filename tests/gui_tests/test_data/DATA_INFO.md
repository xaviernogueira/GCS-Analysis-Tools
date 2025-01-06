# Test Data Info

For the tests here, we are looking at a subset of "Upper North Fork Matilija Creek" (ID=17585756).

**Test User Inputs:** A end-user will start most analyses with only NAIP imagery, and a LAZ file(s). 
Note all test data is stored in `\tests\gui_tests\test_data\`.
* `m_3411930_sw_11_060_20220513_C.tif`: Is a clipped version of four band NAIP imagery from [EarthExplorer](https://earthexplorer.usgs.gov/).
* `USGS_LPC_CA_SoCAL_Wildfires_2018_D18_w2106n1533.laz` (and a second w/ `w2107`): Raw LAZ data covering a piece of the creek.
    * See metadata here: https://www.sciencebase.gov/catalog/item/5db1c935e4b0b0c58b575d88
    * The CRS is encoded in XML medatadata, previously I would download shapefile metdata, but that server seems to be down.
    * Therefore we can use the following manually generated shapefile as the `spatial_shp` argument in `lidar_prep()`.
