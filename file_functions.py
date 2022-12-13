import os
from typing import List, Tuple, Union
from tkinter import *
import subprocess
import logging
import arcpy
import csv

arcpy.env.overwriteOutput = True


def init_logger(filename) -> None:
    """Initializes logger"""
    logging.basicConfig(
        filename=os.path.basename(filename).replace(
            '.py',
            '.log',
        ),
        filemode='w',
        level=logging.INFO,
    )
    stderr_logger = logging.StreamHandler()
    stderr_logger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logging.getLogger().addHandler(stderr_logger)
    return


def cmd(command) -> None:
    """Executes command prompt command"""
    try:
        res = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except Exception:
        msg = 'Command failed: %s' % command
        logging.error(msg)
        raise Exception(msg)

    msg = res.communicate()[1]
    msg_str = str(msg, 'utf-8')

    # if using for LAStools, get rid of the annoying LAStools licensing message.
    # 'Please note that LAStools is not "free" (see http://lastools.org/LICENSE.txt) contact martin.isenburg@rapidlasso.com to clarify licensing terms if needed.',

    if 'http://lastools.org/LICENSE.txt' not in msg_str and len(msg_str) > 0:
        logging.info(msg)
    return


def err_info(func) -> function:
    """Wrapper to show error message when a command fails"""
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            logging.error(e)
    return wrapper


def spatial_license(func) -> function:
    def wrapper(*args, **kwargs):
        arcpy.CheckOutExtension('Spatial')
        func(*args, **kwargs)
        arcpy.CheckInExtension('Spatial')
    return wrapper


def check_use(filepath) -> None:
    """Checks if a file or list of files is in use by another process
    If the file cannot be opened or there is an associated .lock file, it throws an exception.
    """

    if type(filepath) == list:
        for f in filepath:
            check_use(f)
        return

    file_object = None
    if os.path.exists(filepath):
        error_msg = f'{filepath} is open in another program. Close the file and try again.'
        try:
            buffer_size = 8
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'a', buffer_size)
            if file_object:
                for filename in os.listdir(os.path.dirname(filepath)):
                    if filename.startswith(os.path.basename(filepath)) and filename.endswith('.lock'):
                        raise PermissionError(error_msg)

        except IOError:
            raise PermissionError(error_msg)

        finally:
            if file_object:
                file_object.close()
    return


def table_to_csv(
    input_table,
    csv_filepath,
    fld_to_remove_override: List[str] = [],
    keep_fields: List[str] = [],
) -> str:
    """Returns the file path of a csv containing the attributes table of a shapefile or other table"""
    fld_list = arcpy.ListFields(input_table)
    fld_names = [str(fld.name) for fld in fld_list]

    # Either delete specified fields, or only keep specified fields
    if len(fld_to_remove_override) > 0:
        for field in fld_to_remove_override:
            try:
                fld_names.remove(field)
            except Exception:
                logging.warning("Can't delete field: %s" % field)

    elif len(keep_fields) > 0:
        fld_names = [i for i in fld_names if i in keep_fields]

    with open(csv_filepath, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fld_names)
        with arcpy.da.SearchCursor(input_table, fld_names) as cursor:
            for row in cursor:
                writer.writerow(row)
        logging.info(csv_filepath + " CREATED")
    csv_file.close()

    return csv_filepath


def delete_gis_files(file_loc) -> None:
    """This function accepts a GIS file location (eg. \\shapefile.shp) and deletes the file as well
    as any other related file (eg. shapefile.prj, shapefile.cpg). This function supports .tif, .shp, and .dbf"""
    suffix = file_loc[-4:]
    prefix = file_loc[:-4]
    if suffix == '.shp':
        suf_list = [suffix, '.cpg', '.dbf', '.prj',
                    '.sbn', '.sbx', '.shp.xlm', '.shx']

    elif suffix == '.tif':
        suf_list = [suffix, '.tif.aux.xml', '.tfw',
                    '.tif.ovr', '.tif.vat.cpg', '.tif.vat.dbf']

    elif suffix == '.dbf':
        suf_list = [suffix, '.cpg', '.dbf.xml']

    elif suffix == '.csv':
        suf_list = [suffix]
    else:
        suf_list = []

    counter = 0
    for suf in suf_list:
        file = prefix + suf
        if os.path.exists(file):
            try:
                os.remove(file)
            except Exception:
                logging.warning("Couldn't delete %s" % file)
        else:
            counter += 1

    logging.warning(
        f'Couldnt find {counter} files sub-files. '
        'Not normally and issue but if overwrite errors raise this could be the culprit!'
    )


def find_suffix(csv_location) -> str:
    """This function takes a csv table location and finds the suffix unaffected by stage.
    Ex: C://documents//2p3ft_gcs_table.csv would return ft_gcs_table as a string"""
    base = os.path.basename(csv_location)

    if str.isnumeric(base[0]):
        index = 0
        base_snip = base[0]
        while base_snip != 'f' and base_snip != 'm':
            index += 1
            base_snip = base[index]

        suffix = str(base[index:])

    else:
        raise FileNotFoundError(
            'csv filename not suitable. Please have stage height and units in name at the start of the filename. '
            'Ex: 2p3ft_gcs_table.csv or 1m_gcs_table.csv'
        )
    return suffix


def float_keyz_format(z) -> str:
    """This function takes a float key z argument and returns its equivalent formatted string.
    ex: 5.3 -> 5p3, or 10.0 -> 10p0"""

    z_str = ''
    if isinstance(z, float):
        z1, z2 = str(z).split('.', 1)
        z_str = f'{z1}p{z2}'
    elif isinstance(z, int):
        z_str = str(z) + 'p0'

    try:
        return z_str
    except z_str == '':
        raise ValueError(
            'Key z list parameters not valid. Please fill list with int or float.')


def string_to_list(
    string: str,
    format: str = '',
) -> List[Union[None, str]]:
    """Splits a string at each comma and produces a list. format parameter allows the type of the output to
    be designated"""
    out_list = []
    if string != '':
        str_split = string.split(',')
        if format == 'int':
            out_list = [int(i) for i in str_split]
        elif format == 'float':
            out_list = [float(i) for i in str_split]
    return out_list


def get_label_units(projected_file) -> Tuple[str, str, arcpy.SpatialReference]:
    """Input: Any file with a projected spatial reference
    Returns: list storing unit label for filenames, linear unit itself, and the spatial reference file
    Example return: ['m', 'Meters', *spatial_ref_object*]"""

    spatial_ref = arcpy.Describe(projected_file).spatialReference
    unit = spatial_ref.linearUnitName

    if unit == 'Meter':
        u = 'm'
    else:
        u = 'ft'

    return (u, unit, spatial_ref)
