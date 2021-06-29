import arcpy
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')
import file_functions
from file_functions import *
import logging


def least_cost_centerline(DEM, source):
    '''returns a rough centerline using least cost path from source'''
    check_use([DEM, source])

    # make directory for output files
    outdir = os.path.dirname(DEM) + '\\centerline\\'
    if os.path.isdir(outdir) == False:
        os.mkdir(outdir)
        logging.info('Created output directory %s' % outdir)

    ###############################
    # fill sinks in DEM raster
    logging.info('Filling sinks in DEM...')
    check_use(outdir + 'filled_DEM.tif')
    filled_DEM = arcpy.sa.Fill(DEM.split('.aux')[0])
    filled_DEM.save(outdir + 'filled_DEM.tif')
    logging.info('OK')

    ###############################
    # make flow direction raster
    logging.info('Computing flow direction...')
    check_use(outdir + 'flow_dir.tif')
    flow_dir = arcpy.sa.FlowDirection(filled_DEM)
    flow_dir.save(outdir + 'flow_dir.tif')
    logging.info('OK')

    ###############################
    # create least cost path
    logging.info('Computing least cost path...')
    check_use(outdir + 'lc_path.tif')
    lc_path_raster = arcpy.sa.CostPath(source,
                                       filled_DEM,
                                       flow_dir,
                                       path_type='BEST_SINGLE',
                                       destination_field='Id'
                                       )
    lc_path_raster.save(outdir + 'lc_path.tif')
    logging.info('OK')

    ###############################
    # create rough center polyline
    logging.info('Creating rough centerline...')
    check_use(outdir + 'rough_centerline.shp')
    rough_centerline = arcpy.RasterToPolyline_conversion(lc_path_raster,
                                                         outdir + 'rough_centerline.shp',
                                                         simplify='NO_SIMPLIFY'
                                                         )
    logging.info('OK')

    logging.info('Deleting intermediate files...')
    del_files = [filled_DEM, flow_dir, lc_path_raster]
    for f in del_files:
        arcpy.Delete_management(f)
    logging.info('OK')

    return rough_centerline.getOutput(0)


def remove_spurs(line, spur_length=2):
    '''Removes spurs from line smaller than spur_length'''
    check_use([line, line.replace('.shp', '_rm_spurs.shp')])

    logging.info('Removing spurs smaller than % s units...' % str(spur_length))
    # measure lengths
    arcpy.AddGeometryAttributes_management(line, 'LENGTH')
    # convert to layer
    lyr = arcpy.MakeFeatureLayer_management(line, line.replace('.shp', '.lyr'))
    # select the line segments with length > spur_length
    no_spurs = arcpy.SelectLayerByAttribute_management(lyr, where_clause='LENGTH > %s' % str(spur_length))
    # copy this selection to a new feature
    no_spurs = arcpy.CopyFeatures_management(no_spurs, line.replace('.shp', '_rm_spurs.lyr'))
    # delete input line and layers
    del_files = [line, line.replace('.shp', '.lyr'), line.replace('.shp', '_rm_spurs.lyr')]
    for f in del_files:
        arcpy.Delete_management(f)

    logging.info('OK.')

    return line.replace('.shp', '_rm_spurs.shp')


def smooth_centerline(rough_centerline, smooth_distance):
    '''returns a smoothed version of the rough centerline'''

    outdir = os.path.dirname(rough_centerline) + '/'
    check_use([rough_centerline, outdir + 'smooth_centerline.shp'])

    # dissolve rough centerline in case multiple pieces were made...
    centerline = arcpy.Dissolve_management(rough_centerline, outdir + 'dis_rough_centerline.shp')

    logging.info('Smoothing centerline...')
    centerline = arcpy.cartography.SmoothLine(centerline,
                                              outdir + 'smooth_centerline.shp',
                                              algorithm='PAEK',
                                              tolerance=smooth_distance
                                              )

    arcpy.Delete_management(outdir + 'dis_rough_centerline.shp')
    logging.info('OK')

    return centerline.getOutput(0)


def clip_centerline(centerline, channel):
    '''removes parts of centerline falling outside channel polygon'''

    outdir = os.path.dirname(centerline) + '/'
    check_use([centerline, channel, outdir + 'centerline.shp'])

    logging.info('Clipping centerline to channel...')
    centerline = arcpy.Clip_analysis(centerline,
                                     channel,
                                     outdir + 'centerline.shp'
                                     )
    logging.info('OK')
    return centerline.getOutput(0)


@err_info
@spatial_license
def make_centerline(DEM, channel, source, smooth_distance):
    '''Main function for creating, smoothing, and clipping a centerline

    Args:
        DEM: DEM raster
        channel: river channel polygon (for clipping centerline)
        source: flow source polygon (with Id attribute set to 1)
        smooth_distance: centerline smoothness param (tolerance for PAEK smoothing algorithm)

    Returns:
        centerline shapefile
    '''
    outdir = os.path.dirname(DEM) + '/centerline/'
    check_use([DEM,
               channel,
               source,
               outdir + 'filled_DEM.tif',
               outdir + 'flow_dir.tif',
               outdir + 'lc_path.tif',
               outdir + 'rough_centerline.shp',
               outdir + 'smooth_centerline.shp',
               outdir + 'clipped_centerline.shp'
               ])
    rough_centerline = least_cost_centerline(DEM, source)
    rough_centerline = remove_spurs(rough_centerline)
    centerline = smooth_centerline(rough_centerline, smooth_distance)
    centerline = clip_centerline(centerline, channel)

    logging.info('Deleting intermediate files...')
    del_files = [rough_centerline, outdir + 'smooth_centerline.shp']
    for f in del_files:
        arcpy.Delete_management(f)
    logging.info('OK')

    logging.info('Finished: %s' % centerline)
    return centerline


if __name__ == '__main__':
    '''
    #test data
    #input DEM
    DEM = 'G:\Kenny_Rainbow_Basin\DEM\Rainbow_Basin_DEM2.tif'
    #channel polygon
    channel = 'G:\Kenny_Rainbow_Basin\Channel1.shp'
    #flow source polygon
    source = 'G:\Kenny_Rainbow_Basin\DEM\centerline\lc_source.shp'
    #centerline smoothness
    smooth_distance = 20
    '''

    init_logger(__file__)

    # make the GUI window
    root = Tk()
    root.wm_title('Create centerline')

    # specify relevant directories/files

    L1 = Label(root, text='DEM:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='file', ftypes=[('Raster', '*.tif'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='Channel Polygon:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E2.grid(row=1, column=2)
    b2 = Button(root, text='Browse', command=lambda: browse(root, E2, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b2.grid(sticky=W, row=1, column=3)

    L3 = Label(root, text='Flow Source Polygon:')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E3.grid(row=2, column=2)
    b3 = Button(root, text='Browse', command=lambda: browse(root, E3, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b3.grid(sticky=W, row=2, column=3)

    L4 = Label(root, text='Smoothing Distance:')
    L4.grid(sticky=E, row=8, column=1)
    E4 = Entry(root, bd=5)
    E4.insert(END, 20)
    E4.grid(row=8, column=2)

    b = Button(root, text='   Run    ', command=lambda: make_centerline(DEM=E1.get(),
                                                                        channel=E2.get(),
                                                                        source=E3.get(),
                                                                        smooth_distance=float(E4.get())
                                                                        )
               )
    b.grid(sticky=W, row=9, column=2)
    root.grid_rowconfigure(9, minsize=80)

    root.mainloop()

arcpy.CheckInExtension('Spatial')
