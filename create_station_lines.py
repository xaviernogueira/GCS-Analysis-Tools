import arcpy
from arcpy import env
import os
import tkinter
from tkinter import *
from file_functions import *
import logging
arcpy.env.overwriteOutput = True

@err_info
def create_station_lines_function(line_shp, spacing, xs_length, stage=[]):
    '''
    Creates station lines perpendicular to line_shp with given longitudinal spacing and lateral XS length
    (lengths are in units of input line coordinate system)
    '''
    
    # convert line to route feature
    # *** include coordinate priority so we automatically have stationing oriented DOWNSTREAM
    if isinstance(stage, int):
        stage = [stage]
    elif isinstance(stage, float):
        stage = [int(stage)]
    logging.info('Converting input line to route...')
    
    line_fields = [field.name for field in arcpy.ListFields(line_shp)]
    id = [x for x in ['Id', 'arcid', 'ObjectID'] if x in line_fields][0]
    if id==[]:
        raise Exception('Couldn\'t find Id Field in %s' % line_shp)

    route = arcpy.CreateRoutes_lr(line_shp, id, line_shp.replace('.shp', '_rt.shp'))
        
    logging.info('OK.')

    route_id, total_length = 0, 0
    # get route id and total length of route (only considering case of one line for input shapefile)
    for row in arcpy.da.SearchCursor(route, [id, 'SHAPE@LENGTH']):
        route_id = row[0]
        total_length = row[1]

    # total number of XS's
    num_xs = int(total_length / spacing)

    # make route events tables
    logging.info('Creating route events table...')
    locations = [spacing*x for x in range(num_xs)]*2
    offsets = [xs_length*1.0/2]*num_xs+[-xs_length*1.0/2]*num_xs
    route_ids = [route_id]*2*num_xs

    df = pd.DataFrame.from_dict(dict(zip(['LOCATION', 'OFFSET', 'ROUTE'], [locations, offsets, route_ids])))
    event_table = line_shp.replace('.shp', '_ret.csv')
    df.to_csv(event_table, index=False)
    logging.info('OK.')

    # make route event layer
    logging.info('Making route event layer...')
    el_u = arcpy.MakeRouteEventLayer_lr(route, id, event_table, 'ROUTE POINT LOCATION', 'event_layer', 'OFFSET')
    # merge w/ itself to index
    el_name = line_shp.replace('.shp', '_el.shp')
    el = arcpy.Merge_management(el_u, el_name)
    logging.info('OK.')

    # convert to output polyline
    logging.info('Converting points to polyline output...')
    if len(stage) > 0:
        out_name = line_shp.replace('.shp', ('_XS_%sft.shp' % (int(spacing))))
    else:
        out_name = line_shp.replace('.shp', ('_XS_%sx%sft.shp' % (int(spacing), int(xs_length))))
    arcpy.PointsToLine_management(el, out_name, 'LOCATION')
    logging.info('OK.')

    # delete intermediate files
    delfiles = [route, event_table, el_u, el]
    for delfile in delfiles:
        arcpy.Delete_management(delfile)

    logging.info('Finished. Output: %s' % out_name)
    if len(stage) > 0:
        print('Station lines for stage %sft made. Output: %s' % (stage[0], out_name))
    return out_name


if __name__ == '__main__':

    init_logger(__file__)

    # make the GUI
    root = Tk()
    root.wm_title('Create station lines')

    # specify relevant directories/files

    L1 = Label(root, text='centerline:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse',
                command=lambda: browse(root, E1, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                        ('All files', '*')]
                                       )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='spacing:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, 3)
    E2.grid(row=1, column=2)

    L3 = Label(root, text='XS length:')
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5)
    E3.insert(END, 100)
    E3.grid(row=2, column=2)

    L4 = Label(root, text='Stage')
    L4.grid(sticky=E, row=3, column=1)
    E4 = Entry(root, bd=5)
    E4.insert(END, 10)
    E4.grid(row=3, column=2)

    b = Button(root, text='   Run    ',
               command=lambda: create_station_lines_function(line_shp=E1.get(), spacing=float(E2.get()), xs_length=float(E3.get()), stage=[int(E4.get())])
               )
    b.grid(sticky=W, row=4, column=2)
    root.grid_rowconfigure(4, minsize=80)

    root.mainloop()
