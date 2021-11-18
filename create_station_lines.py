import arcpy
from arcpy import env
import os
import tkinter
from tkinter import *
from file_functions import *
import logging

arcpy.env.overwriteOutput = True

@err_info
def create_station_lines_function(line_shp, spacing, xs_length):
    """Creates station lines perpendicular to line_shp with given longitudinal spacing and lateral XS length
    (lengths are in units of input line coordinate system)."""

    # convert line to route feature
    # *** include coordinate priority so we automatically have stationing oriented DOWNSTREAM

    line_dir = os.path.dirname(line_shp)
    init_logger(__file__)  # Initiate log file

    # Initiate temp files folder
    temp_files = os.path.dirname(line_dir) + '\\temp_files'

    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    line_fields = [field.name for field in arcpy.ListFields(line_shp)]
    id = [x for x in ['Id', 'arcid', 'ObjectID'] if x in line_fields][0]
    if id == []:
        raise Exception(r'Couldn\'t find Id Field in %s' % line_shp)

    route = arcpy.CreateRoutes_lr(line_shp, id, line_shp.replace('.shp', '_rt.shp'))
    route_id, total_length = 0, 0

    # get route id and total length of route (only considering case of one line for input shapefile)
    for row in arcpy.da.SearchCursor(route, [id, 'SHAPE@LENGTH']):
        route_id = row[0]
        total_length = row[1]

    # total number of XS's
    num_xs = int(total_length / spacing)

    # make route events tables
    locations = [spacing * x for x in range(num_xs)] * 2
    offsets = [xs_length * 1.0 / 2] * num_xs + [-xs_length * 1.0 / 2] * num_xs
    route_ids = [route_id] * 2 * num_xs

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

    line_name = os.path.basename(line_shp)
    out_name = temp_files + '\\%s' % line_name.replace('.shp', '_XS.shp')

    arcpy.PointsToLine_management(el, out_name, 'LOCATION')
    logging.info('OK.')

    # delete intermediate files
    delfiles = [route, event_table, el_u, el]
    for delfile in delfiles:
        try:
            arcpy.Delete_management(delfile)
        except:
            print('')

    logging.info('Finished. Output: %s' % out_name)

    return out_name


