import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import itertools
import tkinter
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
import subprocess
import logging
import arcpy
import csv

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')
logger = logging.getLogger(__name__)


def init_logger(filename):
    """Initializes logger"""
    logging.basicConfig(filename=os.path.basename(filename).replace('.py', '.log'), filemode='w', level=logging.INFO)
    stderr_logger = logging.StreamHandler()
    stderr_logger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logging.getLogger().addHandler(stderr_logger)
    return


def cmd(command):
    """Executes command prompt command"""
    try:
        res = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        msg = 'Command failed: %s' % command
        logger.error(msg)
        raise Exception(msg)

    msg = res.communicate()[1]

    # if using for LAStools, get rid of the annoying LAStools licensing message.
    if r'http://lastools.org/LICENSE.txt' in msg:
        msg = ''
    # msg = msg.replace(
    # 'Please note that LAStools is not "free" (see http://lastools.org/LICENSE.txt) contact martin.isenburg@rapidlasso.com to clarify licensing terms if needed.',
    # '')
    logger.info(msg)
    return


def browse(root, entry, select='file', ftypes=[('All files', '*')]):
    """GUI button command: opens browser window and adds selected file/folder to entry"""
    if select == 'file':
        filename = filedialog.askopenfilename(parent=root, title='Choose a file', filetypes=ftypes)
        if filename != None:
            entry.delete(0, END)
            entry.insert(END, filename)

    elif select == 'files':
        files = filedialog.askopenfilenames(parent=root, title='Choose files', filetypes=ftypes)
        l = root.tk.splitlist(files)
        entry.delete(0, END)
        entry.insert(END, l)

    elif select == 'folder':
        dirname = filedialog.askdirectory(parent=root, initialdir=entry.get(), title='Choose a directory')
        if len(dirname) > 0:
            entry.delete(0, END)
            entry.insert(END, dirname + '/')


def err_info(func):
    """Wrapper to show error message when a command fails"""
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            logger.info(e)
            # showerror('Error', e)

    return wrapper


def spatial_license(func):
    def wrapper(*args, **kwargs):
        arcpy.CheckOutExtension('Spatial')
        func(*args, **kwargs)
        arcpy.CheckInExtension('Spatial')

    return wrapper


def check_use(filepath):
    """Checks if a file or list of files is in use by another process
    If the file cannot be opened or there is an associated .lock file, it throws an exception.
    """

    if type(filepath) == list:
        for f in filepath:
            check_use(f)
        return

    file_object = None
    if os.path.exists(filepath):
        try:
            buffer_size = 8
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'a', buffer_size)
            if file_object:
                for filename in os.listdir(os.path.dirname(filepath)):
                    if filename.startswith(os.path.basename(filepath)) and filename.endswith('.lock'):
                        logger.error('%s is open in another program. Close the file and try again.' % filepath)
                        raise Exception('%s is open in another program. Close the file and try again.' % filepath)

        except IOError:
            logger.error('%s is open in another program. Close the file and try again.' % filepath)
            raise Exception('%s is open in another program. Close the file and try again.' % filepath)

        finally:
            if file_object:
                file_object.close()
    return


def split_list(l, break_pts):
    """returns list l split up into sublists at break point indices"""
    l_0 = len(l)
    sl = []

    # Return a list containing the input list if no breakpoints indices selected
    if len(break_pts) == 0:
        return [l]

    # Else splits the list and return a list of sub lists. ADJUST SO IT'S NOT BP INDICES BUT RATHER LOCATION VALUES?
    else:
        for brk in break_pts:
            delta_l = l_0 - len(l)
            sl.append(l[:brk - delta_l])
            l = l[brk - delta_l:]
        sl.append(l)
    return sl


def split_reaches(l, new_reach_pts):
    """splits l into sections where new_reach_pts contains the starting indices for each slice"""
    new_reach_pts = sorted(new_reach_pts)
    sl = [l[i1:i2] for i1, i2 in zip(new_reach_pts, new_reach_pts[1:])]
    last_index = new_reach_pts[-1]
    sl.append(l[last_index:])
    return sl


def tif_to_poly(tif):
    """Converts .tif raster to a single polygon covering area that is not null"""
    ras = arcpy.Raster(tif)

    # Make all non-null cells have a value of 1
    int_raster = arcpy.sa.Con(arcpy.sa.IsNull(ras) == False, 1)

    # Covert raster to polygon
    poly = arcpy.RasterToPolygon_conversion(int_raster, tif.replace('.tif', '.shp'), 'NO_SIMPLIFY')
    return poly.getOutput(0)


def tableToCSV(input_table, csv_filepath, fld_to_remove_override=[], keep_fields=[]):
    """Returns the file path of a csv containing the attributes table of a shapefile or other table"""
    fld_list = arcpy.ListFields(input_table)
    fld_names = [str(fld.name) for fld in fld_list]

    # Either delete specified fields, or only keep specified fields
    if len(fld_to_remove_override) > 0:
        for field in fld_to_remove_override:
            try:
                fld_names.remove(field)
            except:
                "Can't delete field: %s" % field

    elif len(keep_fields) > 0:
        fld_names = [i for i in fld_names if i in keep_fields]

    with open(csv_filepath, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fld_names)
        with arcpy.da.SearchCursor(input_table, fld_names) as cursor:
            for row in cursor:
                writer.writerow(row)
        print(csv_filepath + " CREATED")
    csv_file.close()

    return csv_filepath


def delete_gis_files(file_loc):
    """This function accepts a GIS file location (eg. \\shapefile.shp) and deletes the file as well
    as any other related file (eg. shapefile.prj, shapefile.cpg). This function supports .tif, .shp, and .dbf"""
    suffix = file_loc[-4:]
    prefix = file_loc[:-4]
    if suffix == '.shp':
        suf_list = ['.shp', '.cpg', '.dbf', '.prj', '.sbn', '.sbx', '.shp.xlm', '.shx']

    elif suffix == '.tif':
        suf_list = ['.tif', '.tif.aux.xml', '.tfw', '.tif.ovr', '.tif.vat.cpg', '.tif.vat.dbf']

    elif suffix == '.dbf':
        suf_list = ['.dbf', '.cpg', '.dbf.xml']

    elif suffix == '.csv':
        suf_list = ['.csv']

    counter = 0
    for suf in suf_list:
        file = prefix + suf
        if os.path.exists(file):
            try:
                os.remove(file)
            except:
                print("Couldn't delete %s" % file)
        else:
            counter += 1

    print(
        'Couldnt find %s files sub-files. Not normally and issue but if overwrite errors raise this could be the culprit!' % counter)


def find_suffix(csv_location):
    """This function takes a csv table location and finds the suffix unaffected by stage.
    Ex: C://documents//2p3ft_gcs_table.csv would return ft_gcs_table as a string"""
    base = os.path.basename(csv_location)

    if str.isnumeric(base[0]) == True:
        index = 0
        base_snip = base[0]
        while base_snip != 'f' and base_snip != 'm':
            index += 1
            base_snip = base[index]

        suffix = str(base[index:])

    else:
        print(
            'csv filename not suitable. Please have stage height and units in name at the start of the filename. Ex: 2p3ft_gcs_table.csv or 1m_gcs_table.csv')

    return suffix


def float_keyz_format(z):
    '''This function takes a float key z argument and retrusn its equivalent formatted string.
    ex: 5.3 -> 5p3, or 10.0 -> 10p0'''

    z_str = ''
    if z >= 10.0 and isinstance(z, float):
        z_str = (str(z)[0:2] + 'p' + str(z)[3])
    elif z < 10.0 and isinstance(z, float):
        z_str = (str(z)[0] + 'p' + str(z)[2])
    elif isinstance(z, int):
        z_str = str(z) + 'p0'

    try:
        return z_str
    except z_str == '':
        print('Key z list parameters not valid. Please fill list with int or float.')
