import logging
from file_functions import *
import arcpy

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')


def least_cost_centerline(
    DEM: str,
    source: str,
) -> str:
    """returns a rough centerline using least cost path from source"""
    check_use([DEM, source])

    try:
        # make directory for output files
        outdir = os.path.dirname(DEM) + '\\'

        # fill sinks in DEM raster
        logging.info('Filling sinks in DEM...')
        check_use(outdir + 'filled_DEM.tif')
        filled_DEM = arcpy.sa.Fill(DEM.split('.aux')[0])
        filled_DEM.save(outdir + 'filled_DEM.tif')
        logging.info('OK')

        # make flow direction raster
        logging.info('Computing flow direction...')
        check_use(outdir + 'flow_dir.tif')
        flow_dir = arcpy.sa.FlowDirection(filled_DEM)
        flow_dir.save(outdir + 'flow_dir.tif')
        logging.info('OK')

        # create least cost path
        logging.info('Computing least cost path...')
        check_use(outdir + 'lc_path.tif')
        lc_path_raster = arcpy.sa.CostPath(
            source,
            filled_DEM,
            flow_dir,
            path_type='BEST_SINGLE',
            destination_field='Id',
        )
        lc_path_raster.save(outdir + 'lc_path.tif')
        logging.info('OK')

        ###############################
        # create rough center polyline
        logging.info('Creating rough centerline...')
        check_use(outdir + 'rough_centerline.shp')
        rough_centerline = arcpy.RasterToPolyline_conversion(
            outdir + 'lc_path.tif',
            outdir + 'rough_centerline.shp',
            simplify='NO_SIMPLIFY',
        )
        logging.info('OK')

        logging.info('Deleting intermediate files...')
        del_files = [filled_DEM, flow_dir, lc_path_raster]

        for f in del_files:
            try:
                arcpy.Delete_management(f)
            except Exception:
                logging.warning(f'Could not delete intermediate file {f}')
        logging.info('OK')

    except arcpy.ExecuteError:
        logging.info(str(arcpy.GetMessages()))

    return rough_centerline.getOutput(0)


def remove_spurs(
    line: str,
    spur_length: Union[float, int] = 10,
) -> str:
    """Removes spurs from line smaller than spur_length"""
    check_use([line, line.replace('.shp', '_rm_spurs.shp')])

    try:
        logging.info('Removing spurs smaller than % s units...' %
                     str(spur_length))

        # measure lengths
        arcpy.AddGeometryAttributes_management(line, 'LENGTH')

        # convert to layer
        lyr = arcpy.MakeFeatureLayer_management(
            line,
            line.replace('.shp', '.lyr'),
        )

        # select the line segments with length > spur_length
        no_spurs = arcpy.SelectLayerByAttribute_management(
            lyr,
            where_clause=f'LENGTH > {str(spur_length)}',
        )

        # copy this selection to a new feature
        no_spurs = arcpy.CopyFeatures_management(
            no_spurs,
            line.replace('.shp', '_rm_spurs.lyr'),
        )

        # delete input line and layers
        del_files = [
            line.replace('.shp', '.lyr'),
            line.replace('.shp', '_rm_spurs.lyr'),
        ]

        for f in del_files:
            try:
                arcpy.Delete_management(f)
            except Exception:
                logging.warning(f'Could not delete intermediate file {f}')

        logging.info('OK.')

    except arcpy.ExecuteError:
        logging.info(str(arcpy.GetMessages()))

    return line.replace('.shp', '_rm_spurs.shp')


def smooth_centerline(
    rough_centerline: str,
    smooth_distance: Union[float, int],
) -> str:
    """Returns a smoothed version of the rough centerline"""

    outdir = os.path.dirname(rough_centerline) + '\\'
    check_use([rough_centerline, outdir + 'no_clip_thalweg_centerline.shp'])

    # dissolve rough centerline in case multiple pieces were made...
    centerline = arcpy.Dissolve_management(
        rough_centerline,
        outdir + 'dis_rough_centerline.shp',
    )

    logging.info('Smoothing centerline...')
    centerline = arcpy.cartography.SmoothLine(
        centerline,
        outdir + 'no_clip_thalweg_centerline.shp',
        algorithm='PAEK',
        tolerance=smooth_distance,
    )

    try:
        del_file = outdir + 'dis_rough_centerline.shp'
        arcpy.Delete_management(del_file)
    except Exception:
        logging.warning(f'Could not delete intermediate file {del_file}')
    del del_file
    logging.info('OK')

    return centerline.getOutput(0)


def clip_centerline(
    centerline: str,
    clip_poly1: str,
    clip_poly2: str = '',
) -> str:
    """Removes parts of centerline falling outside channel polygon"""
    outdir = os.path.dirname(centerline) + '\\'
    check_use([centerline, clip_poly1, outdir + 'thalweg_centerline.shp'])

    if clip_poly2 == '':
        logging.info('Clipping centerline to channel...')
        centerline = arcpy.Clip_analysis(
            centerline,
            clip_poly1,
            outdir + 'thalweg_centerline.shp',
        )
        logging.info('OK')
    else:
        logging.info('Clipping centerline to channel and lidar_extent...')
        centerline = arcpy.Clip_analysis(
            centerline,
            clip_poly1,
            outdir + 'thalweg_centerline1.shp',
        )

        centerline = arcpy.Clip_analysis(
            centerline,
            clip_poly2,
            outdir + 'thalweg_centerline.shp',
        )

        del_file = outdir + 'thalweg_centerline1.shp'
        try:
            arcpy.Delete_management(del_file)
        except Exception:
            logging.warning(f'Could not delete intermediate file {del_file}')
        logging.info('OK')

    return centerline.getOutput(0)


@err_info
@spatial_license
def make_centerline(
    DEM: str,
    channel: str,
    lidar_extent: str,
    source: str,
    smooth_distance: Union[float, int],
) -> str:
    """Main function for creating, smoothing, and clipping a centerline

    Args:
        DEM: DEM raster
        channel: river channel polygon (for clipping centerline)
        source: flow source polygon (with Id attribute set to 1)
        smooth_distance: centerline smoothness param (tolerance for PAEK smoothing algorithm)

    Returns:
        centerline shapefile
    """
    outdir = os.path.dirname(DEM)

    # arcpy.env.extent = lidar_extent - causes a sneaky exception but work around may need to be introduced

    check_use(
        [
            DEM,
            channel,
            source,
            lidar_extent,
            outdir + '\\filled_DEM.tif',
            outdir + '\\flow_dir.tif',
            outdir + '\\lc_path.tif',
            outdir + '\\rough_centerline.shp',
            outdir + '\\smooth_centerline.shp',
            outdir + '\\thalweg_centerline.shp',
        ],
    )

    rough_centerline = least_cost_centerline(DEM, source)
    rough_centerline = remove_spurs(rough_centerline)
    centerline = smooth_centerline(rough_centerline, smooth_distance)
    centerline = clip_centerline(centerline, channel, lidar_extent)

    logging.info('Deleting intermediate files...')
    del_files = [
        rough_centerline, outdir + '\\rough_centerline.shp',
        outdir + '\\no_clip_thalweg_centerline.shp',
    ]

    for f in del_files:
        try:
            arcpy.Delete_management(f)
        except Exception:
            logging.warning(f'Could not delete intermediate file {f}')

    logging.info('Finished: %s' % centerline)
    return centerline
