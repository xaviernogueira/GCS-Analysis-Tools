from tkinter import *
from file_functions import *
import os
import sys
import shutil
import numpy as np
import logging
import gc

gc.collect()


##########################
# first let's define some functions that will be helpful

def las_files(directory):
    """returns list of all .las/.laz files in directory (at top level)"""
    l = []
    for name in os.listdir(directory):
        if name.endswith('.las') or name.endswith('.laz'):
            l.append(directory + name)

    return l


# input working directory for LAStools and directory containing .las/.laz files
# creates a .txt file for LAStools containing list of .las/.laz file names
# returns the name of the .txt file.
def lof_text(pwd, src):
    """creates a .txt file in pwd (LAStools bin) containing a list of .las/.laz filenames from src directory"""
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


# input .las/.laz filename, outputs point density (after running lasinfo)
def pd(filename):
    """returns point density from lasinfo output .txt file"""
    # name of txt output file from lasinfo
    filename = filename[:-4] + '.txt'
    f = open(filename, 'r')
    text = f.readlines()
    for line in text:
        if line.startswith('point density:'):
            d = line.split(' ')
            d = d[d.index('returns') + 1]
            return float(d)


def get_largest(directory):
    """returns name of largest file in directory"""
    largest_so_far = 0
    filename = ''
    for name in os.listdir(directory):
        size = os.path.getsize(os.path.join(directory, name))
        if size > largest_so_far:
            largest_so_far = size
            filename = name

    return os.path.join(directory, filename)


def pts(filename, lastoolsdir):
    """returns number of points in las file"""

    # call lasinfo on the file
    cmd('%slasinfo.exe -i %s -otxt -histo number_of_returns 1' % (lastoolsdir, filename))
    # name of txt output file from lasinfo
    txt = filename[:-4] + '.txt'
    f = open(txt, 'r')
    text = f.readlines()
    for line in text:
        if line.endswith('element(s)\n'):
            d = line.split(' ')
            d = d[d.index('for') + 1]
            return int(d)


# the main function that runs when 'run' button is clicked
@err_info
def process_lidar(lastoolsdir,
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
                  fine_offset
                  ):
    """Executes main LAStools processing workflow. See readme for more info."""

    # initialize the logger
    init_logger(__file__)

    classes = ['01-Default',
               '02-Ground',
               '05-Vegetation',
               '06-Building'
               ]

    if (ground_poly != '') and (keep_orig_pts == True):
        # run on coarse and fine settings, need to clip and remove duplicates after merging
        outdirs = ['00_separated',
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
                   '13_veg_rm_duplicates'
                   ]

    elif (ground_poly == '') and (keep_orig_pts == True):
        # only classify with coarse settings, no clipping, but need to remove duplicates
        outdirs = ['00_separated',
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
                   '13_veg_rm_duplicates'
                   ]

    elif (ground_poly == '') and (keep_orig_pts == False):
        # only classify with coarse setting, no clipping or removing duplicates necessary
        outdirs = ['00_separated',
                   '00_declassified',
                   '01_tiled',
                   '02_lasground_new',
                   '03_lasheight',
                   '04_lasclassify',
                   '05_lastile_rm_buffer',
                   '06_separated'
                   ]

    elif (ground_poly != '') and (keep_orig_pts == False):
        # run on coarse and fine settings, clip, but no removing duplicates needed
        outdirs = ['00_separated',
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
                   '12_veg_merged'
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
        sepdirs = [lidardir + '00_separated',
                   lidardir + '06a_separated_coarse',
                   lidardir + '06b_separated_fine'
                   ]
    else:
        sepdirs = [lidardir + '00_separated',
                   lidardir + '06_separated'
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
                lidar_files.append(path + name)  # Used to have a '/' added between path and name

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
    cmd('%slasinfo.exe -lof %s -set_classification 1 -otxt -cd' % (lastoolsdir, lof))

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
        ds.append(pd(filename))
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
    logging.info('Checking if largest file has < 1.5M pts (to avoid licensing restrictions)...')
    largest_file = get_largest(odir)
    num = pts(largest_file, lastoolsdir)
    if num < 1500000:
        logging.info('Largest file has %i points, tiles small enough.' % num)
    else:
        logging.info('Tile size not small enough. Retrying with a smaller tile size...')
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
            logging.info('Checking if largest file has < 1.5M pts (to avoid licensing restrictions)...')
            largest_file = get_largest(odir)
            num = pts(largest_file, lastoolsdir)
            if num >= 1500000:
                logging.info('Tile size not small enough. Retrying with a smaller tile size...')

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
                odir
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

    cmd('%slasheight.exe -lof %s -cores %i -odir %s -olas' % (lastoolsdir, lof, cores, odir))

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '02b_lasground_new_fine/')
        odir = lidardir + '03b_lasheight_fine/'

        cmd('%slasheight.exe -lof %s -cores %i -odir %s -olas' % (lastoolsdir, lof, cores, odir))

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

    cmd('%slasclassify.exe -lof %s -cores %i %s -odir %s -olas' % (lastoolsdir, lof, cores, units_code, odir))

    logging.info('OK')

    if ground_poly != '':
        logging.info('Classifying non-ground points on fine setting...')

        lof = lof_text(lastoolsdir, lidardir + '03b_lasheight_fine/')
        odir = lidardir + '04b_lasclassify_fine/'

        cmd('%slasclassify.exe -lof %s -cores %i %s -odir %s -olas' % (lastoolsdir, lof, cores, units_code, odir))

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
    cmd('%slastile.exe -lof %s -cores %i -remove_buffer -odir %s -olas' % (lastoolsdir, lof, cores, odir))

    if ground_poly != '':
        lof = lof_text(lastoolsdir, lidardir + '04b_lasclassify_fine/')
        odir = lidardir + '05b_lastile_rm_buffer_fine/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slastile.exe -lof %s -cores %i -remove_buffer -odir %s -olas' % (lastoolsdir, lof, cores, odir))

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
        logging.info('Clipping ground points to inverse ground polygon on coarse setting...')

        # keep points outside ground polygon for coarse setting (-interior flag)
        lof = lof_text(lastoolsdir, lidardir + '06a_separated_coarse' + '/' + '02-Ground' + '/')
        odir = lidardir + '07a_ground_clipped_coarse/'

        cmd('%slasindex.exe -lof %s -cores %i' % (lastoolsdir, lof, cores))
        cmd('%slasclip.exe -lof %s -cores %i -poly %s -interior -donuts -odir %s -olas' % (
            lastoolsdir, lof, cores, ground_poly, odir))

        logging.info('OK')

        logging.info('Clipping ground points to ground polygon on fine setting...')

        # keep points inside ground polygon for fine setting
        lof = lof_text(lastoolsdir, lidardir + '06b_separated_fine' + '/' + '02-Ground' + '/')
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
        sources = [lidardir + '07a_ground_clipped_coarse/', lidardir + '07b_ground_clipped_fine/']

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
        cmd('%slasduplicate.exe -lof %s -cores %i -lowest_z -odir %s -olas' % (lastoolsdir, lof, cores, odir))

        logging.info('OK')

    ##########################
    # merge new veg points

    if ground_poly != '':
        logging.info('Merging new vegetation points from coarse and fine run...')

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
            sources = [lidardir + '11_veg_new_clipped/', lidardir + '00_separated' + '/' + '05-Vegetation' + '/']
        else:
            sources = [lidardir + '06_separated' + '/' + '05-Vegetation' + '/',
                       lidardir + '00_separated' + '/' + '05-Vegetation' + '/']
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
        cmd('%slasduplicate.exe -lof %s -cores %i -lowest_z -odir %s -olas' % (lastoolsdir, lof, cores, odir))

        logging.info('OK')

    logging.info('Processing finished.')
    logging.info('Outputs in:')
    logging.info(ground_results)
    logging.info('')
    logging.info(veg_results)

    return


#####################################################################

if __name__ == '__main__':
    # initialize the logger
    init_logger(__file__)

    # make the GUI window
    root = Tk()
    root.wm_title('LiDAR Reprocessing App (based on LAStools)')

    # specify relevant directories/files

    L1 = Label(root, text='LAStools /bin/ directory:')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse', command=lambda: browse(root, E1, select='folder'))
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='LiDAR data directory:')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.insert(END, '/'.join(sys.path[0].split('\\')[:-1]) + '/')
    E2.grid(row=1, column=2)
    b2 = Button(root, text='Browse', command=lambda: browse(root, E2, select='folder'))
    b2.grid(sticky=W, row=1, column=3)

    L3 = Label(root, text='Ground area .shp file (optional):')
    shp_var = StringVar()
    L3.grid(sticky=E, row=2, column=1)
    E3 = Entry(root, bd=5, textvariable=shp_var)
    E3.grid(row=2, column=2)
    b3 = Button(root, text='Browse', command=lambda: browse(root, E3, select='file', ftypes=[('Shapefile', '*.shp'),
                                                                                             ('All files', '*')]
                                                            )
                )
    b3.grid(sticky=W, row=2, column=3)


    # if no ground shapefile is provided, disable the fine setting and just run on "coarse"
    def trace_choice(*args):
        if shp_var.get() == '':
            for widget in [E1b, E2b, E3b, E4b, E5b]:
                widget.config(state=DISABLED)
        else:
            for widget in [E1b, E2b, E3b, E4b, E5b]:
                widget.config(state='normal')


    shp_var.trace('w', trace_choice)

    # specify lasground_new parameters

    root.grid_rowconfigure(5, minsize=80)

    LC1 = Label(root, text='standard/coarse classification parameters:')
    LC1.grid(row=5, column=0, columnspan=2)

    L1a = Label(root, text='step size:')
    L1a.grid(sticky=E, row=6)
    E1a = Entry(root, bd=5)
    E1a.grid(row=6, column=1)

    L2a = Label(root, text='bulge:')
    L2a.grid(sticky=E, row=7)
    E2a = Entry(root, bd=5)
    E2a.grid(row=7, column=1)

    L3a = Label(root, text='spike:')
    L3a.grid(sticky=E, row=8)
    E3a = Entry(root, bd=5)
    E3a.grid(row=8, column=1)

    L4a = Label(root, text='down spike:')
    L4a.grid(sticky=E, row=9)
    E4a = Entry(root, bd=5)
    E4a.grid(row=9, column=1)

    L5a = Label(root, text='offset:')
    L5a.grid(sticky=E, row=10)
    E5a = Entry(root, bd=5)
    E5a.grid(row=10, column=1)

    LC2 = Label(root, text='fine classification parameters (in ground area):')
    LC2.grid(row=5, column=2, columnspan=2)

    L1b = Label(root, text='step size:')
    L1b.grid(sticky=E, row=6, column=2)
    E1b = Entry(root, bd=5, state=DISABLED)
    E1b.grid(row=6, column=3)

    L2b = Label(root, text='bulge:')
    L2b.grid(sticky=E, row=7, column=2)
    E2b = Entry(root, bd=5, state=DISABLED)
    E2b.grid(row=7, column=3)

    L3b = Label(root, text='spike:')
    L3b.grid(sticky=E, row=8, column=2)
    E3b = Entry(root, bd=5, state=DISABLED)
    E3b.grid(row=8, column=3)

    L4b = Label(root, text='down spike:')
    L4b.grid(sticky=E, row=9, column=2)
    E4b = Entry(root, bd=5, state=DISABLED)
    E4b.grid(row=9, column=3)

    L5b = Label(root, text='offset:')
    L5b.grid(sticky=E, row=10, column=2)
    E5b = Entry(root, bd=5, state=DISABLED)
    E5b.grid(row=10, column=3)

    # specify units
    L5 = Label(root, text='Units')
    L5.grid(sticky=W, row=11, column=2)
    root.grid_rowconfigure(11, minsize=80)
    unit_var = StringVar()
    R5m = Radiobutton(root, text='Meters', variable=unit_var, value=' ')
    R5m.grid(sticky=E, row=12, column=1)
    R5f = Radiobutton(root, text='US Feet', variable=unit_var, value=' -feet -elevation_feet ')
    R5f.grid(row=12, column=2)
    unit_var.set(' ')

    # specify number of cores
    L4 = Label(root, text='Number of cores for processing')
    L4.grid(sticky=E, row=13, column=1, columnspan=2)
    root.grid_rowconfigure(13, minsize=80)
    core_num = IntVar()
    R1 = Radiobutton(root, text='1', variable=core_num, value=1)
    R1.grid(sticky=E, row=14, column=1)
    R2 = Radiobutton(root, text='2', variable=core_num, value=2)
    R2.grid(row=14, column=2)
    R4 = Radiobutton(root, text='4', variable=core_num, value=4)
    R4.grid(sticky=W, row=14, column=3)
    R8 = Radiobutton(root, text='8', variable=core_num, value=8)
    R8.grid(sticky=E, row=15, column=1)
    R16 = Radiobutton(root, text='16', variable=core_num, value=16)
    R16.grid(row=15, column=2)
    R32 = Radiobutton(root, text='32', variable=core_num, value=32)
    R32.grid(sticky=W, row=15, column=3)
    core_num.set(16)

    L5 = Label(root, text='Keep original ground/veg points: ')
    L5.grid(sticky=E, row=16, column=1)
    keep_originals = BooleanVar()
    C1 = Checkbutton(root, variable=keep_originals)
    C1.grid(sticky=W, row=16, column=2)
    keep_originals.set(True)

    # make 'Run' button in GUI to call the process_lidar() function
    b = Button(root, text='    Run    ', command=lambda: process_lidar(lastoolsdir=E1.get(),
                                                                       lidardir=E2.get(),
                                                                       ground_poly=E3.get(),
                                                                       cores=core_num.get(),
                                                                       units_code=unit_var.get()[1:-1],
                                                                       keep_orig_pts=keep_originals.get(),
                                                                       coarse_step=E1a.get(),
                                                                       coarse_bulge=E2a.get(),
                                                                       coarse_spike=E3a.get(),
                                                                       coarse_down_spike=E4a.get(),
                                                                       coarse_offset=E5a.get(),
                                                                       fine_step=E1b.get(),
                                                                       fine_bulge=E2b.get(),
                                                                       fine_spike=E3b.get(),
                                                                       fine_down_spike=E4b.get(),
                                                                       fine_offset=E5b.get()
                                                                       )
               )

    b.grid(sticky=W, row=17, column=2)
    root.grid_rowconfigure(17, minsize=80)

    root.mainloop()
