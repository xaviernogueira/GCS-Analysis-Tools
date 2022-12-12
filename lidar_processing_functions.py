from tkinter import *
from file_functions import *
import os
import shutil
import numpy as np
import logging
import gc

gc.collect()


##########################
# first let's define some functions that will be helpful

def las_files(
    directory: str,
) -> List[str]:
    """returns list of all .las/.laz files in directory (at top level)"""
    l = []
    for name in os.listdir(directory):
        if name.endswith('.las') or name.endswith('.laz'):
            l.append(directory + name)

    return l


def lof_text(
    pwd: str,
    src: str,
) -> str:
    """Creates a .txt file in pwd (LAStools bin) containing a list of .las/.laz filenames from src directory

    Args:
        pwd: LAStools bin path
        src: working directory path for LAStools and directory containing .las/.laz files

    Returns:
        The name of the .txt file.
    """
    filename = pwd + 'file_list.txt'

    f = open(filename, 'w+')

    if type(src) == str:
        for i in las_files(src):
            f.write('%s\n' % i)
    else:
        # this is the case when there are multiple source folders
        for i in [name for source in src for name in las_files(source)]:
            f.write('%s\n' % i)

    f.close()
    return filename


def point_density(
    filename: str,
) -> float:
    """Returns point density from lasinfo output .txt file via the .las/.laz filename"""

    # name of txt output file from lasinfo
    filename = filename[:-4] + '.txt'
    f = open(filename, 'r')
    text = f.readlines()
    for line in text:
        if line.startswith('point density:'):
            d = line.split(' ')
            d = d[d.index('returns') + 1]
            return float(d)


def get_largest(
    directory: str,
) -> str:
    """returns name of largest file in directory"""
    largest_so_far = 0
    filename = ''
    for name in os.listdir(directory):
        size = os.path.getsize(os.path.join(directory, name))
        if size > largest_so_far:
            largest_so_far = size
            filename = name

    return os.path.join(directory, filename)


def pts(
    filename: str,
    lastoolsdir: str,
) -> int:
    """returns number of points in las file"""

    # call lasinfo on the file
    cmd('%slasinfo.exe -i %s -otxt -histo number_of_returns 1' %
        (lastoolsdir, filename))
    # name of txt output file from lasinfo
    txt = filename[:-4] + '.txt'
    f = open(txt, 'r')
    text = f.readlines()
    for line in text:
        if line.endswith('element(s)\n'):
            d = line.split(' ')
            d = d[d.index('for') + 1]
            return int(d)


@err_info
def process_lidar(
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
) -> None:
    """Executes main LAStools processing workflow. See readme for more info."""

    classes = [
        '01-Default',
        '02-Ground',
        '05-Vegetation',
        '06-Building',
    ]

    if (ground_poly != '') and (keep_orig_pts == True):
        # run on coarse and fine settings, need to clip and remove duplicates after merging
        outdirs = [
            '00_separated',
            '00_declassified',
            '01_tiled',
            '02a_lasground_new_coarse',
            '02b_lasground_new_fine',
            '03a_lasheight_coarse',
            '03b_lasheight_fine',
            '04a_lasclassify_coarse',
            '04b_lasclassify_fine',
            '05a_lastile_rm_buffer_coarse',
            '05b_lastile_rm_buffer_fine',
            '06a_separated_coarse',
            '06b_separated_fine',
            '07a_ground_clipped_coarse',
            '07b_ground_clipped_fine',
            '08_ground_merged',
            '09_ground_rm_duplicates',
            '10_veg_new_merged',
            '11_veg_new_clipped',
            '12_veg_merged',
            '13_veg_rm_duplicates',
        ]

    elif (ground_poly == '') and (keep_orig_pts == True):
        # only classify with coarse settings, no clipping, but need to remove duplicates
        outdirs = [
            '00_separated',
            '00_declassified',
            '01_tiled',
            '02_lasground_new',
            '03_lasheight',
            '04_lasclassify',
            '05_lastile_rm_buffer',
            '06_separated',
            '08_ground_merged',
            '09_ground_rm_duplicates',
            '12_veg_merged',
            '13_veg_rm_duplicates',
        ]

    elif (ground_poly == '') and (keep_orig_pts == False):
        # only classify with coarse setting, no clipping or removing duplicates necessary
        outdirs = [
            '00_separated',
            '00_declassified',
            '01_tiled',
            '02_lasground_new',
            '03_lasheight',
            '04_lasclassify',
            '05_lastile_rm_buffer',
            '06_separated',
        ]

    elif (ground_poly != '') and (keep_orig_pts == False):
        # run on coarse and fine settings, clip, but no removing duplicates needed
        outdirs = [
            '00_separated',
            '00_declassified',
            '01_tiled',
            '02a_lasground_new_coarse',
            '02b_lasground_new_fine',
            '03a_lasheight_coarse',
            '03b_lasheight_fine',
            '04a_lasclassify_coarse',
            '04b_lasclassify_fine',
            '05a_lastile_rm_buffer_coarse',
            '05b_lastile_rm_buffer_fine',
            '06a_separated_coarse',
            '06b_separated_fine',
            '07a_ground_clipped_coarse',
            '07b_ground_clipped_fine',
            '08_ground_merged',
            '10_veg_new_merged',
            '11_veg_new_clipped',
            '12_veg_merged',
        ]

    # make new directories for output from each step in processing
    for ind, outdir in enumerate(outdirs):
        if os.path.isdir(lidardir + outdir) == False:
            os.mkdir(lidardir + outdir)

    if len(os.listdir(lidardir + outdirs[0])) != 0:
        msg = 'Output directories must initially be empty. Move or delete the data currently in output directories.'
        logging.error(msg)
        raise Exception(msg)

    # in each 'separated' folder, create subdirs for each class type
    if ground_poly != '':
        sepdirs = [
            lidardir + '00_separated',
            lidardir + '06a_separated_coarse',
            lidardir + '06b_separated_fine',
        ]
    else:
        sepdirs = [
            lidardir + '00_separated',
            lidardir + '06_separated',
        ]
    for sepdir in sepdirs:
        for class_type in classes:
            class_dir = sepdir + '/' + class_type
            if os.path.isdir(class_dir) == False:
                os.mkdir(class_dir)

    logging.info('Created directories for output data')

    ##################################
    # create declassified points

    logging.info('Declassifying copy of original point cloud...')

    # get list of filenames for original LiDAR data (all .las and .laz files in lidardir)
    lidar_files = []
    for path, subdirs, files in os.walk(lidardir):
        for name in files:
            if name.endswith('.las') or name.endswith('.laz'):
                # Used to have a '/' added between path and name
                lidar_files.append(path + name)

    if lidar_files == []:
        msg = 'No .las or .laz files in %s or its subdirectories' % lidardir
        logging.error(msg)
        raise Exception(msg)

    # remove duplicate files already copied to declassified in re-run scenarios
    for name in lidar_files:
        if '00_declassified' in name:
            lidar_files.remove(name)

    # copy original files into '00_declassified' folder
    out_fol = lidardir + '00_declassified/'
    for name in lidar_files:
        name_base = os.path.basename(name)

        if name_base not in os.listdir(out_fol):
            shutil.copyfile(name, out_fol + name_base)

    # make list of files for LASTools to process, which is output
    lof = lof_text(lastoolsdir, lidardir + '00_declassified/')

    # call LAStools command to declassify points and get point density
    cmd('%slasinfo.exe -lof %s -set_classification 1 -otxt -cd' %
        (lastoolsdir, lof))

    logging.info('OK')

    # separate original data by class type

    logging.info('Separating original data by class type...')

    filename = lastoolsdir + 'file_list.txt'
    f = open(filename, 'w+')
    for i in lidar_files:
        f.write('%s\n' % i)
    f.close()
    lof = filename

    for class_type in classes:
        odir = lidardir + '00_separated' + '/' + class_type + '/'
        class_code = int(class_type.split('-')[0])
        cmd('%slas2las.exe -lof %s -cores %i -keep_classification %i -odir %s -olas' % (
            lastoolsdir, lof, cores, class_code, odir))

    logging.info('OK')

    ########################
    # create tiling (max 1.5M pts per tile)

    logging.info('Creating tiling...')

    # get point density for each .las file
    ds = []
    for filename in las_files(lidardir + '00_declassified/'):
        ds.append(point_density(filename))
    # use max point density out of all files to determine tile size
    max_d = max(ds)

    # width of square tile so we have max of 1.5M pts per tile (assuming same number of points per tile)
    # throw in another factor of 0.5 to make sure tiles will be small enough, round to nearest 10
    tile_size = round(0.5 * np.sqrt((1.5 * 10 ** 6) / max_d), -1)

    logging.info('Using tile size of %i' % tile_size)

    odir = (lidardir + '01_tiled/')

    # call LAStools command to create tiling
    cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
    cmd('%slastile.exe -lof %s -cores %i -tile_size %i -buffer 5 -faf -odir %s -o tile.las -olas' % (
        lastoolsdir, lof, cores, tile_size, odir))

    logging.info('OK')

    # check to make sure tiles are small enough
    logging.info(
        'Checking if largest file has < 1.5M pts (to avoid licensing restrictions)...')
    largest_file = get_largest(odir)
    num = pts(largest_file, lastoolsdir)
    if num < 1500000:
        logging.info('Largest file has %i points, tiles small enough.' % num)
    else:
        logging.info(
            'Tile size not small enough. Retrying with a smaller tile size...')
        while num >= 1500000:
            # delete original set of tiles
            folder = odir
            for the_file in os.listdir(folder):
                file_path = os.path.join(folder, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except:
                    logging.warning('Couldn\'nt delete %s' % file_path)
            # redo tiling
            tile_size = int(tile_size * num * 1.0 / 1500000)
            logging.info('Using tile size of %i' % tile_size)

            cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
            cmd('%slastile.exe -lof %s -cores %i -o tile.las -tile_size %i -buffer 5 -faf -odir %s -olas' % (
                lastoolsdir, lof, cores, tile_size, odir))
            # recheck largest tile number of points
            logging.info(
                'Checking if largest file has < 1.5M pts (to avoid licensing restrictions)...')
            largest_file = get_largest(odir)
            num = pts(largest_file, lastoolsdir)
            if num >= 1500000:
                logging.info(
                    'Tile size not small enough. Retrying with a smaller tile size...')

    logging.info('OK')

    ########################
    # run lasground_new on coarse and fine settings

    logging.info('Running ground classification on coarse setting...')

    lof = lof_text(lastoolsdir, lidardir + '01_tiled/')

    if ground_poly != '':
        odir = lidardir + '02a_lasground_new_coarse/'
    else:
        odir = lidardir + '02_lasground_new'

    cmd(
        '%slasground_new.exe -lof %s -cores %i %s -step %s -bulge %s -spike %s -down_spike %s -offset %s -hyper_fine -odir %s -olas' % (
            lastoolsdir,
            lof,
            cores,
            units_code,
            coarse_step,
            coarse_bulge,
            coarse_spike,
            coarse_down_spike,
            coarse_offset,
            odir
        )
    )

    logging.info('OK')

    if ground_poly != '':
        logging.info('Running ground classification on fine setting...')

        odir = lidardir + '02b_lasground_new_fine/'

        cmd(
            '%slasground_new.exe -lof %s -cores %i %s -step %s -bulge %s -spike %s -down_spike %s -offset %s -hyper_fine -odir %s -olas' % (
                lastoolsdir,
                lof,
                cores,
                units_code,
                fine_step,
                fine_bulge,
                fine_spike,
                fine_down_spike,
                fine_offset,
                odir,
            )
        )

        logging.info('OK')

    ##########################
    # run lasheight on each data set

    logging.info('Measuring height above ground for non-ground points...')

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '02a_lasground_new_coarse/')
        odir = lidardir + '03a_lasheight_coarse/'
    else:
        lof = lof_text(lastoolsdir, lidardir + '02_lasground_new/')
        odir = lidardir + '03_lasheight/'

    cmd('%slasheight.exe -lof %s -cores %i -odir %s -olas' %
        (lastoolsdir, lof, cores, odir))

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '02b_lasground_new_fine/')
        odir = lidardir + '03b_lasheight_fine/'

        cmd('%slasheight.exe -lof %s -cores %i -odir %s -olas' %
            (lastoolsdir, lof, cores, odir))

    logging.info('OK')

    ##########################
    # run lasclassify on each data set

    logging.info('Classifying non-ground points on coarse setting...')

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '03a_lasheight_coarse/')
        odir = lidardir + '04a_lasclassify_coarse/'
    else:
        lof = lof_text(lastoolsdir, lidardir + '03_lasheight/')
        odir = lidardir + '04_lasclassify/'

    cmd('%slasclassify.exe -lof %s -cores %i %s -odir %s -olas' %
        (lastoolsdir, lof, cores, units_code, odir))

    logging.info('OK')

    if ground_poly != '':
        logging.info('Classifying non-ground points on fine setting...')

        lof = lof_text(lastoolsdir, lidardir + '03b_lasheight_fine/')
        odir = lidardir + '04b_lasclassify_fine/'

        cmd('%slasclassify.exe -lof %s -cores %i %s -odir %s -olas' %
            (lastoolsdir, lof, cores, units_code, odir))

        logging.info('OK')

    ##########################
    # remove tile buffers on each data set

    logging.info('Removing tile buffers...')

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '04a_lasclassify_coarse/')
        odir = lidardir + '05a_lastile_rm_buffer_coarse/'
    else:
        lof = lof_text(lastoolsdir, lidardir + '04_lasclassify/')
        odir = lidardir + '05_lastile_rm_buffer/'

    cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
    cmd('%slastile.exe -lof %s -cores %i -remove_buffer -odir %s -olas' %
        (lastoolsdir, lof, cores, odir))

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '04b_lasclassify_fine/')
        odir = lidardir + '05b_lastile_rm_buffer_fine/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slastile.exe -lof %s -cores %i -remove_buffer -odir %s -olas' %
            (lastoolsdir, lof, cores, odir))

    logging.info('OK')

    ##########################
    # separate into files for each class type

    logging.info('Separating points by class type on coarse setting...')

    # coarse
    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '05a_lastile_rm_buffer_coarse/')
        podir = lidardir + '06a_separated_coarse'
    else:
        lof = lof_text(lastoolsdir, lidardir + '05_lastile_rm_buffer/')
        podir = lidardir + '06_separated'

    for class_type in classes:
        odir = podir + '/' + class_type + '/'
        class_code = int(class_type.split('-')[0])

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slas2las.exe -lof %s -cores %i -keep_classification %i -odir %s -olas' % (
            lastoolsdir, lof, cores, class_code, odir))

    logging.info('OK')

    if (ground_poly) == '' and (keep_orig_pts == False):
        ground_results = podir + '/' + '02-Ground' + '/'
        veg_results = podir + '/' + '05-Vegetation' + '/'

    if ground_poly != '':
        logging.info('Separating points by class type on fine setting...')

        # fine
        lof = lof_text(lastoolsdir, lidardir + '05b_lastile_rm_buffer_fine/')

        for class_type in classes:
            odir = lidardir + '06b_separated_fine' + '/' + class_type + '/'
            class_code = int(class_type.split('-')[0])

            cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
            cmd('%slas2las.exe -lof %s -cores %i -keep_classification %i -odir %s -olas' % (
                lastoolsdir, lof, cores, class_code, odir))

        logging.info('OK')

    ##########################
    # clip ground data sets with ground polygon
    if ground_poly != '':
        logging.info(
            'Clipping ground points to inverse ground polygon on coarse setting...')

        # keep points outside ground polygon for coarse setting (-interior flag)
        lof = lof_text(lastoolsdir, lidardir +
                       '06a_separated_coarse' + '/' + '02-Ground' + '/')
        odir = lidardir + '07a_ground_clipped_coarse/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasclip.exe -lof %s -cores %i -poly %s -interior -donuts -odir %s -olas' % (
            lastoolsdir, lof, cores, ground_poly, odir))

        logging.info('OK')

        logging.info(
            'Clipping ground points to ground polygon on fine setting...')

        # keep points inside ground polygon for fine setting
        lof = lof_text(lastoolsdir, lidardir +
                       '06b_separated_fine' + '/' + '02-Ground' + '/')
        odir = lidardir + '07b_ground_clipped_fine/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasclip.exe -lof %s -cores %i -poly %s -donuts -odir %s -olas' % (
            lastoolsdir, lof, cores, ground_poly, odir))

        logging.info('OK')

    ##########################
    # merge

    # merge processed ground points with original data set ground points

    if keep_orig_pts:
        logging.info('Merging new and original ground points...')
        if ground_poly != '':
            sources = [lidardir + '07a_ground_clipped_coarse/', lidardir + '07b_ground_clipped_fine/',
                       lidardir + '00_separated' + '/' + '02-Ground' + '/']
        else:
            sources = [lidardir + '06_separated' + '/' + '02-Ground' + '/',
                       lidardir + '00_separated' + '/' + '02-Ground' + '/']
    # just use new points
    elif ground_poly != '':
        logging.info('Merging new ground points...')
        sources = [lidardir + '07a_ground_clipped_coarse/',
                   lidardir + '07b_ground_clipped_fine/']

    if (keep_orig_pts == True) or (ground_poly != ''):
        lof = lof_text(lastoolsdir, sources)
        odir = lidardir + '08_ground_merged/'
        ground_results = odir  # will be overwritten if rm_duplicates block runs

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slastile.exe -lof %s -cores %i -o tile.las -tile_size %i -faf -odir %s -olas' % (
            lastoolsdir, lof, cores, tile_size, odir))

        logging.info('OK')

    ##########################
    # remove duplicate ground points

    if keep_orig_pts:
        logging.info('Removing duplicate ground points...')
        lof = lof_text(lastoolsdir, lidardir + '08_ground_merged/')
        odir = lidardir + '09_ground_rm_duplicates/'
        ground_results = odir

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasduplicate.exe -lof %s -cores %i -lowest_z -odir %s -olas' %
            (lastoolsdir, lof, cores, odir))

        logging.info('OK')

    ##########################
    # merge new veg points

    if ground_poly != '':
        logging.info(
            'Merging new vegetation points from coarse and fine run...')

        sources = [lidardir + '06a_separated_coarse' + '/' + '05-Vegetation' + '/',
                   lidardir + '06b_separated_fine' + '/' + '05-Vegetation' + '/']
        lof = lof_text(lastoolsdir, sources)
        odir = lidardir + '10_veg_new_merged/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slastile.exe -lof %s -cores %i -o tile.las -tile_size %i -faf -odir %s -olas' % (
            lastoolsdir, lof, cores, tile_size, odir))

        logging.info('OK')

        #########################
        # clip new veg points
        # keeping points outside the ground polygon

        logging.info('Clipping new vegetation points...')

        lof = lof_text(lastoolsdir, lidardir + '10_veg_new_merged/')
        odir = lidardir + '11_veg_new_clipped/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasclip.exe -lof %s -cores %i -poly %s -interior -donuts -odir %s -olas' % (
            lastoolsdir, lof, cores, ground_poly, odir))

        logging.info('OK')

    #########################
    # merge with original veg points

    if (keep_orig_pts == True):
        logging.info('Merging new and original vegetation points...')
        if ground_poly != '':
            sources = [
                lidardir + '11_veg_new_clipped/',
                lidardir + '00_separated' + '/' + '05-Vegetation' + '/',
            ]
        else:
            sources = [
                lidardir + '06_separated' + '/' + '05-Vegetation' + '/',
                lidardir + '00_separated' + '/' + '05-Vegetation' + '/',
            ]
    elif ground_poly != '':
        logging.info('Retiling new vegetation points...')
        sources = [lidardir + '11_veg_new_clipped/']

    if (keep_orig_pts == True) or (ground_poly != ''):
        lof = lof_text(lastoolsdir, sources)
        odir = lidardir + '12_veg_merged/'
        veg_results = odir  # will be overwritten if rm_duplicates block runs

    cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
    cmd('%slastile.exe -lof %s -cores %i -o tile.las -tile_size %i -faf -odir %s -olas' % (
        lastoolsdir, lof, cores, tile_size, odir))

    logging.info('OK')

    #########################
    # remove duplicate veg points
    if keep_orig_pts:
        logging.info('Removing duplicate vegetation points...')

        lof = lof_text(lastoolsdir, lidardir + '12_veg_merged/')
        odir = lidardir + '13_veg_rm_duplicates/'
        veg_results = odir
        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasduplicate.exe -lof %s -cores %i -lowest_z -odir %s -olas' %
            (lastoolsdir, lof, cores, odir))

        logging.info('OK')

    logging.info('Processing finished.')
    logging.info('Outputs in:')
    logging.info(ground_results)
    logging.info('')
    logging.info(veg_results)

    return
