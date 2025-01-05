
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
    # First item defines XS length and spacing, second item described smoothing distance
    params = [m_spacing, smooth_dist]

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
