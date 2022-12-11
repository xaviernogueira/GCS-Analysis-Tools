import os
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk
from lidar_to_DEM_functions import lidar_footprint, define_ground_polygon, \
    lidar_to_raster, detrend_prep
from lidar_processing_functions import process_lidar
from dem_detrending_functions import fit_params_txt, make_residual_plot, linear_fit_plot, \
    diagnostic_quick_plot, detrend_that_raster, linear_fit, prep_xl_file
from wetted_area_functions import prep_small_inc, pdf_cdf_plotting, stage_centerlines
from gcs_analysis_functions import extract_gcs
from file_functions import string_to_list
from arcpy import HillShade_3d


class GCSGraphicUserInterface(tk.Frame):

    # Initialize the class, this allows the GUI to run when the code is ran
    def __init__(
            self,
            master=None,
    ):

        # Initialize the tk frame that will hold all tabs holding each processing step
        tk.Frame.__init__(
            self,
            master,
        )
        self.pack()

        # initialize attribute files CHANGE
        self.dem = StringVar()
        self.det_dem = StringVar()
        self.centerline = StringVar()
        self.station_lines = StringVar()

        # set header
        self.master.title('Geomorphic Covariance Structure (GCS) analysis GUI')
        self.master.iconbitmap('win_icon.ico')

        # Create color scheme and padding for the frame/self
        self.bg_color = 'LightSkyBlue1'
        self.padding = 5

        # controls window width
        ww = self.master.winfo_screenwidth() / 2

        # controls window height
        wh = self.master.winfo_screenwidth() / 2

        # position relative to screen width and ww
        wx = (self.master.winfo_screenwidth() - ww) / 2

        # position relative to screen height and wh
        wy = (self.master.winfo_screenheight() - wh) / 2

        # set window height and location
        self.master.geometry("%dx%d+%d+%d" % (ww, wh, wx, wy))

        # set widget styles
        self.style = ttk.Style()

        # Adding the Breeze ttk theme https://github.com/MaxPerl/ttk-Breeze
        breeze_dir = os.getcwd() + '\\ttk-Breeze-master'

        self.tk.call(
            'lappend',
            'auto_path',
            breeze_dir,
        )
        self.tk.call(
            'source',
            breeze_dir + '\\breeze.tcl',
        )
        self.style.theme_use('Breeze')

        # initialize tab handler
        self.tab_container = ttk.Notebook(master)

        self.tab_names = [
            'LiDAR Data prep',
            'DEM generation',
            'Thalweg centerline',
            'Detrend DEM',
            'Flow-stage modeling',
            'GCS analysis',
            'River Builder prep',
        ]

        self.tabs = {}
        for tab_name in self.tab_names:
            tab = ttk.Frame(self.tab_container)
            self.tab_container.add(
                tab,
                text=tab_name,
            )
            self.tab_container.pack(
                expand=1,
                fill="both",
            )
            self.tabs[tab_name] = tab

        # denoted padding between tkinter widgets
        pad = 5

        # Define functions used in multiple windows
        ###################################################################

        def browse(
                root,
                entry,
                select='file',
                ftypes=[('All files', '*')],
        ) -> None:
            """GUI button command: opens browser window and adds selected file/folder to entry"""
            if select == 'file':
                filename = filedialog.askopenfilename(
                    parent=root,
                    title='Choose a file',
                    filetypes=ftypes,
                )
                if filename is not None:
                    entry.delete(0, END)
                    entry.insert(END, filename)

            elif select == 'files':
                files = filedialog.askopenfilenames(
                    parent=root,
                    title='Choose files',
                    filetypes=ftypes,
                )
                l = root.tk.splitlist(files)
                entry.delete(0, END)
                entry.insert(END, l)

            elif select == 'folder':
                dirname = filedialog.askdirectory(
                    parent=root,
                    initialdir=entry.get(),
                    title='Choose a directory',
                )
                if len(dirname) > 0:
                    entry.delete(0, END)
                    # Used to add a \ at the end, may have to bring back if errors occur
                    entry.insert(END, dirname)

        def open_popup(
                title: str,
                image: str,
        ) -> None:
            """Opens a new window showing only an image and a caption displaying image path.
            Inputs: A title that populates the window header, and an image path supported by PIL"""
            self.im = Image.open(image)

            top = Toplevel(root)
            top.geometry()
            top.title(title)

            self.ph = ImageTk.PhotoImage(self.im, master=top)
            self.label = Label(top, image=self.ph)
            self.label.image = self.ph
            self.label.grid(row=1, column=1, columnspan=3)

            self.label2 = Label(top, text='Image saved @ %s' % image)
            self.label2.grid(row=2, column=1)

        # LiDAR prep (filling tabs w/ widgets)
        ######################################################################

        root = self.tabs['LiDAR Data prep']

        def lidar_prep(
                lasbin: str,
                lidardir: str,
                spatial_shp: str,
                naip_folder: str,
                ndvi_thresh: float = 0.4,
                aoi_shp: str = '',
        ) -> None:
            """This function is ran by the LiDAR Data prep tab.
             Outputs: A 1m (or other resolution) DEM modeling bare ground LiDAR returns"""

            print('Unzipping LAZ files...')
            lidar_footprint(
                lasbin,
                lidardir,
                spatial_shp,
            )
            print('Done')

            print('Generating inputs for LiDAR processing...')
            foot = lidardir + '\\las_footprint.shp'
            define_ground_polygon(
                foot,
                lidardir,
                naip_folder,
                ndvi_thresh,
                aoi_shp,
            )
            print('Done')

        self.l_lasbin1 = ttk.Label(
            root,
            text='LAStools /bin/ directory:',
        )
        self.l_lasbin1.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
        )

        self.e_lasbin1 = ttk.Entry(root)
        self.e_lasbin1.insert(END, str(os.getcwd() + '\\LAStools\\bin'))
        self.e_lasbin1.grid(
            row=0,
            column=2,
            pady=pad,
        )

        self.b_lasbin1 = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                entry=self.e_lasbin1,
                select='folder',
            ),
        )
        self.b_lasbin1.grid(
            sticky=W,
            row=0,
            column=3,
            pady=pad,
        )

        self.l_lidardir1 = ttk.Label(
            root,
            text='LiDAR data directory:',
        )
        self.l_lidardir1.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.e_lidardir1 = ttk.Entry(root)
        self.e_lidardir1.insert(END, '')
        self.e_lidardir1.grid(
            row=1,
            column=2,
            pady=pad,
        )
        self.b_lidardir1 = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_lidardir1,
                select='folder',
            ),
        )
        self.b_lidardir1.grid(
            sticky=W,
            row=1,
            column=3,
            pady=pad,
        )

        self.l_in_spatialref = ttk.Label(
            root,
            text='LiDAR spatial reference (.shp):',
        )
        self.l_in_spatialref.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
        )

        self.e_in_spatialref = ttk.Entry(root)
        self.e_in_spatialref.insert(END, '')
        self.e_in_spatialref.grid(
            row=2,
            column=2,
            pady=pad,
        )

        self.b_in_spatialref = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_in_spatialref,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*')
                ],
            ),
        )
        self.b_in_spatialref.grid(
            sticky=W,
            row=2,
            column=3,
            pady=pad,
        )

        self.l_naip = ttk.Label(
            root,
            text='NAIP imagery folder:',
        )
        self.l_naip.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.e_naip = ttk.Entry(root)
        self.e_naip.insert(END, '')
        self.e_naip.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.b_naip = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_naip,
                select='folder',
            ),
        )
        self.b_naip.grid(
            sticky=W,
            row=3,
            column=3,
            pady=pad,
        )

        self.l_ndvi = ttk.Label(
            root,
            text='NDVI vegetation threshold:',
        )
        self.l_ndvi.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )
        self.e_ndvi = ttk.Entry(root)
        self.e_ndvi.grid(
            sticky=E,
            row=4,
            column=2,
            pady=pad,
        )
        self.e_ndvi.insert(0, 0.40)
        self.e_ndvi.grid(
            sticky=E,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_aoi = ttk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_aoi.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
        )

        self.e_aoi = ttk.Entry(root)
        self.e_aoi.insert(END, '')
        self.e_aoi.grid(
            row=5,
            column=2,
            pady=pad,
        )

        self.b_aoi = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_aoi,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_aoi.grid(
            sticky=W,
            row=5,
            column=3,
            pady=pad,
        )

        self.spacer1 = ttk.Label(root, text='')
        self.spacer1.grid(
            sticky=W,
            row=6,
            column=1,
            pady=pad,
        )
        self.prep_run = ttk.Button(
            root,
            text='Run',
            command=lambda: lidar_prep(
                self.e_lasbin1.get(),
                self.e_lidardir1.get(),
                self.e_in_spatialref.get(),
                self.e_naip.get(),
                float(self.e_ndvi.get()),
                self.e_aoi.get(),
            ),
        )
        self.prep_run.grid(
            sticky=E,
            row=6,
            column=2,
        )
        root.grid_rowconfigure(
            17,
            minsize=80,
        )

        self.instruct = ttk.Label(
            root,
            text='     Verify vegetation mask accuracy after running!',
        )
        self.instruct.grid(
            sticky=EW,
            row=6,
            column=3,
            pady=pad,
        )

        input_ref_shp = self.e_in_spatialref.get()

        # DEM generation
        ######################################################################
        # LiDAR processing to DEM generation widgets
        root = self.tabs['DEM generation']

        def dem_generation(
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
                aoi_shp,
                dem_resolution,
                dem_method,
                tri_meth,
                void_meth,
        ):
            """This function used LAStools to generate LAS file tiles prepresenitng bare ground, and then converts them
            to a high resolution DEM (1m resolution is default)"""

            # We carry input spatial ref over from the above process, but we should still convert from shp to ref object
            print('Processing LiDAR to remove vegetation points...')
            las_folder = lidardir + '\\las_files\\'
            process_lidar(
                lastoolsdir + '\\',
                las_folder,
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
            )
            print('Done')

            print('Generating a %sm resolution DEM...' % dem_resolution)
            dem = lidar_to_raster(
                lidardir,
                ground_poly,
                aoi_shp,
                dem_method,
                tri_meth,
                void_meth,
                m_cell_size=float(dem_resolution),
            )
            print('Done')

            print('Generating hillshade raster for the DEM...')
            hill_out = lidardir + '\\hillshade.tif'
            HillShade_3d(
                dem,
                hill_out,
            )
            print('Done')

        self.l_lasbin = ttk.Label(
            root,
            text='LAStools /bin/ directory:',
        )
        self.l_lasbin.grid(
            sticky=E,
            row=0,
            column=1,
        )

        self.e_lasbin = ttk.Entry(root)
        self.e_lasbin.insert(END, os.getcwd() + '\\LAStools\\bin')
        self.e_lasbin.grid(
            row=0,
            column=2,
        )
        self.b_lasbin = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_lasbin,
                select='folder',
            ),
        )
        self.b_lasbin.grid(
            sticky=W,
            row=0,
            column=3,
        )

        self.l_lidardir = ttk.Label(
            root,
            text='LiDAR data directory:',
        )
        self.l_lidardir.grid(
            sticky=E,
            row=1,
            column=1,
        )

        self.e_lidardir = ttk.Entry(root)
        self.e_lidardir.insert(END, '')
        self.e_lidardir.grid(
            row=1,
            column=2,
        )

        self.b_lidardir = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_lidardir,
                select='folder',
            ),
        )
        self.b_lidardir.grid(
            sticky=W,
            row=1,
            column=3,
        )

        self.shp_var = StringVar()

        self.l_ground_shp = ttk.Label(
            root,
            text='Ground polygon (.shp):',
        )
        self.l_ground_shp.grid(
            sticky=E,
            row=2,
            column=1,
        )

        self.e_ground_shp = ttk.Entry(
            root,
            textvariable=self.shp_var,
        )
        self.e_ground_shp.grid(
            row=2,
            column=2,
        )

        self.b_ground_shp = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_ground_shp,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_ground_shp.grid(
            sticky=W,
            row=2,
            column=3,
        )

        self.l_out_spatialref = ttk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_out_spatialref.grid(
            sticky=E,
            row=3,
            column=1,
        )

        self.e_out_spatialref = ttk.Entry(root)
        self.e_out_spatialref.insert(END, '')
        self.e_out_spatialref.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.b_out_spatialref = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_out_spatialref,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_out_spatialref.grid(
            sticky=W,
            row=3,
            column=3,
            pady=pad,
        )

        # if no ground shapefile is provided, disable the fine setting and just run on "coarse"
        def trace_choice(*args):
            fine_entries = [
                self.e_f_step,
                self.e_f_bulge,
                self.e_f_spike,
                self.e_f_dspike,
                self.e_f_offset,
            ]
            if self.shp_var.get() == '':
                for widget in fine_entries:
                    widget.config(state=DISABLED)
            else:
                for widget in fine_entries:
                    widget.config(state='normal')

        self.shp_var.trace('w', trace_choice)

        # specify lasground_new parameters
        root.grid_rowconfigure(5, minsize=80)

        self.l_coarse_class = ttk.Label(
            root,
            text='standard/coarse classification parameters:',
        )
        self.l_coarse_class.grid(
            row=5,
            column=0,
            columnspan=2,
        )

        self.l_c_step = ttk.Label(
            root,
            text='step size:',
        )
        self.l_c_step.grid(
            sticky=E,
            row=6,
        )

        self.e_c_step = ttk.Entry(root)
        self.e_c_step.grid(
            row=6,
            column=1,
        )

        self.l_c_bulge = ttk.Label(
            root,
            text='bulge:',
        )
        self.l_c_bulge.grid(
            sticky=E,
            row=7,
        )

        self.e_c_bulge = ttk.Entry(root)
        self.e_c_bulge.grid(
            row=7,
            column=1,
        )

        self.l_c_spike = ttk.Label(
            root,
            text='spike:',
        )
        self.l_c_spike.grid(
            sticky=E,
            row=8,
        )
        self.e_c_spike = ttk.Entry(root)
        self.e_c_spike.grid(
            row=8,
            column=1,
        )

        self.l_c_dspike = ttk.Label(
            root,
            text='down spike:',
        )
        self.l_c_dspike.grid(
            sticky=E,
            row=9,
        )

        self.e_c_dspike = ttk.Entry(root)
        self.e_c_dspike.grid(
            row=9,
            column=1,
        )

        self.l_c_offset = ttk.Label(
            root,
            text='offset:',
        )
        self.l_c_offset.grid(
            sticky=E,
            row=10,
        )
        self.e_c_offset = ttk.Entry(root)
        self.e_c_offset.grid(
            row=10,
            column=1,
        )

        self.l_fine_class = ttk.Label(
            root,
            text='fine classification parameters (in ground area):',
        )
        self.l_fine_class.grid(
            row=5,
            column=2,
            columnspan=2,
        )

        self.l_f_step = ttk.Label(
            root,
            text='step size:',
        )
        self.l_f_step.grid(
            sticky=E,
            row=6,
            column=2,
        )

        self.e_f_step = ttk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_step.grid(
            row=6,
            column=3,
        )

        self.l_f_bulge = ttk.Label(
            root,
            text='bulge:',
        )
        self.l_f_bulge.grid(
            sticky=E,
            row=7,
            column=2,
        )

        self.e_f_bulge = ttk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_bulge.grid(
            row=7,
            column=3,
        )

        self.l_f_spike = ttk.Label(
            root,
            text='spike:',
        )
        self.l_f_spike.grid(
            sticky=E,
            row=8,
            column=2,
        )

        self.e_f_spike = ttk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_spike.grid(
            row=8,
            column=3,
        )

        self.l_f_dspike = ttk.Label(
            root,
            text='down spike:',
        )
        self.l_f_dspike.grid(
            sticky=E,
            row=9,
            column=2,
        )

        self.e_f_dspike = ttk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_dspike.grid(
            row=9,
            column=3,
        )

        self.l_f_offset = ttk.Label(
            root,
            text='offset:',
        )
        self.l_f_offset.grid(
            sticky=E,
            row=10,
            column=2,
        )

        self.e_f_offset = ttk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_offset.grid(
            row=10,
            column=3,
        )

        # specify units
        self.l_lidar_units = ttk.Label(
            root,
            text='Units:',
        )
        self.l_lidar_units.grid(
            sticky=W,
            row=11,
            column=2,
        )
        root.grid_rowconfigure(11, minsize=30)

        self.lidar_units = StringVar()

        self.r_lidar_meters = ttk.Radiobutton(
            root,
            text='Meters',
            variable=self.lidar_units,
            value=' ',
        )
        self.r_lidar_meters.grid(
            sticky=E,
            row=12,
            column=1,
        )
        self.r_lidar_feet = ttk.Radiobutton(
            root,
            text='US Feet',
            variable=self.lidar_units,
            value=' -feet -elevation_feet ',
        )
        self.r_lidar_feet.grid(
            row=12,
            column=2,
            pady=pad,
        )
        self.lidar_units.set(' ')

        # specify number of cores
        self.l_lidar_cores = ttk.Label(
            root,
            text='Number of cores for processing:',
        )
        self.l_lidar_cores.grid(
            sticky=E,
            row=13,
            column=1,
            columnspan=2,
        )
        root.grid_rowconfigure(13, minsize=30)

        self.core_num = IntVar()

        self.r1_lidar = ttk.Radiobutton(
            root,
            text='1',
            variable=self.core_num, value=1,
        )
        self.r1_lidar.grid(
            sticky=E,
            row=14,
            column=1,
        )

        self.r2_lidar = ttk.Radiobutton(
            root,
            text='2',
            variable=self.core_num,
            value=2,
        )
        self.r2_lidar.grid(
            row=14,
            column=2,
        )

        self.r4_lidar = ttk.Radiobutton(
            root,
            text='4',
            variable=self.core_num,
            value=4,
        )
        self.r4_lidar.grid(
            sticky=W,
            row=14,
            column=3,
        )

        self.r8_lidar = ttk.Radiobutton(
            root,
            text='8',
            variable=self.core_num,
            value=8,
        )
        self.r8_lidar.grid(
            sticky=E,
            row=15,
            column=1,
        )

        self.r16_lidar = ttk.Radiobutton(
            root,
            text='16',
            variable=self.core_num,
            value=16,
        )
        self.r16_lidar.grid(
            row=15,
            column=2,
        )

        self.r32_lidar = ttk.Radiobutton(
            root,
            text='32',
            variable=self.core_num,
            value=32,
        )
        self.r32_lidar.grid(
            sticky=W,
            row=15,
            column=3,
            pady=pad,
        )
        self.core_num.set(16)

        self.l_keep_orig_lidar = ttk.Label(
            root,
            text='Keep original ground/veg points: ',
        )
        self.l_keep_orig_lidar.grid(
            sticky=E,
            row=16,
            column=1,
        )

        self.keep_orig_lidar = BooleanVar()

        self.c_keep_orig_lidar = ttk.Checkbutton(
            root,
            variable=self.keep_orig_lidar,
        )
        self.c_keep_orig_lidar.grid(
            sticky=W,
            row=16,
            column=2,
        )
        self.keep_orig_lidar.set(True)

        units = self.lidar_units.get()[1:-1]

        if units == '':
            default = 1
        else:
            default = 3.28

        self.l_dem_res = ttk.Label(
            root,
            text='DEM resolution (meters):',
        )
        self.l_dem_res.grid(
            sticky=E,
            row=17,
            column=2,
        )

        self.e_dem_res = ttk.Entry(root)
        self.e_dem_res.grid(
            sticky=W,
            row=17,
            column=3,
            pady=pad,
        )
        self.e_dem_res.insert(
            0,
            default,
        )

        # choose binning or triangulation (TIN) based DEM interpolation, be sure to document
        methods = [
            'BINNING',
            'TRIANGULATION',
        ]

        self.l_dem_meth = ttk.Label(
            root,
            text='Select interpolation method:',
        )
        self.l_dem_meth.grid(
            sticky=E,
            row=18,
            column=1,
            pady=pad,
        )

        self.e_dem_meth = StringVar()

        self.option_menu1 = ttk.OptionMenu(
            root,
            self.e_dem_meth,
            *methods,
        )
        self.option_menu1.grid(
            sticky=W,
            row=18,
            column=2,
            pady=pad,
        )

        # select binning method, only relevant if binning is selected as the interpolation
        void_meths = [
            'LINEAR',
            'SIMPLE',
            'NATURAL_NEIGHBOR',
        ]

        self.l_void_meth = ttk.Label(
            root,
            text='Void fill method (for binning:',
        )
        self.l_void_meth.grid(
            sticky=E,
            row=19,
            column=1,
            pady=pad,
        )

        self.e_void_meth = StringVar()

        self.option_menu2 = ttk.OptionMenu(
            root,
            self.e_void_meth,
            *void_meths,
        )
        self.option_menu2.grid(
            sticky=W,
            row=19,
            column=2,
            pady=pad,
        )

        # select binning method, only relevant if binning is selected as the interpolation
        tri_meths = [
            'LINEAR',
            'NATURAL_NEIGHBOR',
        ]

        self.l_tri_meth = ttk.Label(
            root,
            text='Triangulation method:',
        )
        self.l_tri_meth.grid(
            sticky=E,
            row=20,
            column=1,
            pady=pad,
        )
        self.e_tri_meth = StringVar()

        self.option_menu3 = ttk.OptionMenu(
            root,
            self.e_tri_meth,
            *tri_meths,
        )
        self.option_menu3.grid(
            sticky=W,
            row=20,
            column=2,
            pady=pad,
        )

        # make 'Run' ttk.Button in GUI to call the process_lidar() function
        self.b_lidar_run = ttk.Button(
            root,
            text='    Run    ',
            command=lambda: dem_generation(
                lastoolsdir=self.e_lasbin.get(),
                lidardir=self.e_lidardir.get(),
                ground_poly=self.e_ground_shp.get(),
                cores=self.core_num.get(),
                units_code=units,
                keep_orig_pts=self.keep_orig_lidar.get(),
                coarse_step=self.e_c_step.get(),
                coarse_bulge=self.e_c_bulge.get(),
                coarse_spike=self.e_c_spike.get(),
                coarse_down_spike=self.e_c_dspike.get(),
                coarse_offset=self.e_c_offset.get(),
                fine_step=self.e_f_step.get(),
                fine_bulge=self.e_f_bulge.get(),
                fine_spike=self.e_f_spike.get(),
                fine_down_spike=self.e_f_dspike.get(),
                fine_offset=self.e_f_offset.get(),
                aoi_shp=self.e_out_spatialref.get(),
                dem_resolution=self.e_dem_res.get(),
                dem_method=self.e_dem_meth.get(),
                tri_meth=self.e_tri_meth.get(),
                void_meth=self.e_void_meth.get(),
            ),
        )

        self.b_lidar_run.grid(
            sticky=W,
            row=20,
            column=3,
        )
        root.grid_rowconfigure(
            20,
            minsize=40,
        )

        # Generate thalweg centerline and extract elevation profile
        ######################################################################
        root = self.tabs[
            'Thalweg centerline',
        ]

        self.remind = ttk.Label(
            root,
            text='Create upstream flow polygon in ArcMap/Pro',
        )
        self.remind.grid(
            sticky=E,
            row=0,
            column=0,
        )

        self.l_flow_poly = ttk.Label(
            root,
            text='Upstream flow polygon (.shp):',
        )
        self.l_flow_poly.grid(
            sticky=E,
            row=1,
            column=0,
            pady=pad,
        )

        self.e_flow_poly = ttk.Entry(root)
        self.e_flow_poly.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
        )
        self.e_flow_poly.insert(END, '')
        self.e_flow_poly.grid(
            row=1,
            column=1,
            pady=pad,
            padx=5,
        )

        self.b_flow_poly = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_flow_poly,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_flow_poly.grid(
            sticky=W,
            row=1,
            column=2,
            pady=pad,
        )

        self.l_extent = ttk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_extent.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )

        self.e_extent = ttk.Entry(root)
        self.e_extent.grid(
            row=2,
            column=1,
            pady=pad,
        )
        self.e_extent.insert(END, '')
        self.e_extent.grid(
            row=2,
            column=1,
            pady=pad,
            padx=5,
        )

        self.b_extent = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_extent,
                select='file',
                ftypes=[
                    ('Shapefile', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_extent.grid(
            sticky=W,
            row=2,
            column=2,
            pady=pad,
        )

        self.l_filt = ttk.Label(
            root,
            text='Filter passes (15x default):',
        )
        self.l_filt.grid(
            sticky=E,
            row=3,
            column=0,
            pady=pad,
        )

        self.e_filt = ttk.Entry(root)
        self.e_filt.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
        )
        self.e_filt.insert(END, 15)
        self.e_filt.grid(
            row=3,
            column=1,
            pady=pad,
            padx=5,
        )

        self.l_smooth = ttk.Label(
            root,
            text='Smoothing distance (meters):',
        )
        self.l_smooth.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_smooth = ttk.Entry(root)
        self.e_smooth.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )
        self.e_smooth.insert(END, 6)
        self.e_smooth.grid(
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )

        self.l_dem = ttk.Label(
            root,
            text='DEM (.tif):',
        )
        self.l_dem.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_dem = ttk.Entry(root)
        self.e_dem.grid(
            row=5,
            column=1,
            pady=pad,
        )
        self.e_dem.insert(END, '')
        self.e_dem.grid(
            row=5,
            column=1,
            pady=pad,
            padx=5,
        )

        self.b_dem = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_dem,
                select='file',
                ftypes=[
                    ('TIFF', '*.tif'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_dem.grid(
            sticky=W,
            row=5,
            column=2,
            pady=pad,
        )

        # create run botton to create smoothed centerline
        self.b_detrend_prep1 = ttk.Button(
            root,
            text='    Run    ',
            command=lambda: detrend_prep(
                dem=self.e_dem.get(),
                flow_poly=self.e_flow_poly.get(),
                aoi_shp=self.e_extent.get(),
                filt_passes=self.e_filt.get(),
                smooth_dist=self.e_smooth.get(),
                m_spacing=1,
                centerline_verified=False,
            ),
        )
        self.b_detrend_prep1.grid(
            sticky=W,
            row=6,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(6, minsize=50)

        self.l_step = ttk.Label(
            root,
            text='Verify centerline quality (edit if necessary), then run below...',
        )
        self.l_step.grid(
            sticky=E,
            row=7,
            column=0,
        )

        # Create run button to generate a thalweg elevation table from 1m spaced station points
        self.b_detrend_prep2 = ttk.Button(
            root,
            text='    Generate thalweg profile    ',
            command=lambda: detrend_prep(
                dem=self.e_dem.get(),
                flow_poly=self.e_flow_poly.get(),
                aoi_shp=self.e_extent.get(),
                filt_passes=self.e_filt.get(),
                smooth_dist=self.e_smooth.get(),
                m_spacing=1,
                centerline_verified=True,
            ),
        )
        self.b_detrend_prep2.grid(
            sticky=E,
            row=8,
            column=0,
            pady=15,
        )
        root.grid_rowconfigure(
            6,
            minsize=50,
        )

        # DEM detrending
        ######################################################################
        root = self.tabs[
            'Detrend DEM',
        ]

        def show_plot(xyz_csv):
            """Plots thalweg elevation profile is plotted but not saved"""
            # set up directory for plots
            out_dir = os.path.dirname(xyz_csv) + '\\detrending_plots'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            # make arrays. Save and display (pop up) the thalweg elevation plot
            out_list = prep_xl_file(
                xyz_csv,
                in_columns=[
                    'LOCATION',
                    'POINT_X',
                    'POINT_Y',
                    'Value',
                ],
            )

            plot = diagnostic_quick_plot(
                location_np=out_list[0],
                z_np=out_list[1],
                out_dir=out_dir,
            )

            open_popup('Thalweg elevation profile', plot)
            print('Displaying thalweg elevation profile plot!')

        def show_fit_plots(xyz_csv, breakpoints):
            """Linear fit and residuals are plotted but not saved"""
            # set up directory for plots
            out_dir = os.path.dirname(xyz_csv) + '\\detrending_plots'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            # format breakpoints list from string
            breakpoint_list = string_to_list(
                breakpoints,
                format='int',
            )

            # apply linear fit to input csv
            out_list = prep_xl_file(
                xyz_csv,
                in_columns=[
                    'LOCATION',
                    'POINT_X',
                    'POINT_Y',
                    'Value',
                ],
            )

            fit_out = linear_fit(
                location_np=out_list[0],
                z_np=out_list[1],
                xyz_table_loc=xyz_csv,
                bp_list=breakpoint_list,
            )

            # save and display (pop up) linear fit and residual plots. Generate txt
            fit_plot = linear_fit_plot(
                location_np=out_list[0],
                z_np=out_list[1],
                fit_params=fit_out[0],
                fit_np=fit_out[1],
                out_dir=out_dir,
            )

            res_plot = make_residual_plot(
                location_np=out_list[0],
                residual_np=fit_out[2],
                r2=fit_out[3],
                out_dir=out_dir,
            )

            txt = fit_params_txt(
                fit_params=fit_out[0],
                bp_list=breakpoint_list,
                out_dir=out_dir,
            )

            open_popup(
                'Linear fit w/ breakpoints: %s' % breakpoint_list,
                fit_plot,
            )
            open_popup(
                'Residual plot w/ breakpoints: %s' % breakpoint_list,
                res_plot,
            )
            print('Text file listing linear piecewise fit components @ %s' % txt)

        def detrend(
                xyz_csv,
                in_dem,
                aoi_shp
        ) -> None:
            """Detrends the input raster based on .csv stored x, y, z, z_fit values for thalweg points"""

            print('Detrending DEM...')
            detrended_dem = detrend_that_raster(
                xyz_csv=xyz_csv,
                in_dem=in_dem,
                aoi_shp=aoi_shp,
            )
            print('Done')
            print('Detrended DEM @ %s' % detrended_dem)

        self.l_xyz = ttk.Label(
            root,
            text='Thalweg profile csv:',
        )
        self.l_xyz.grid(
            sticky=E,
            row=0,
            column=0,
        )

        self.e_xyz = ttk.Entry(root)
        self.e_xyz.grid(
            stick=E,
            row=0,
            column=1,
        )
        self.e_xyz.insert(END, '')
        self.e_xyz.grid(
            stick=E,
            row=0,
            column=1,
            padx=5,
        )

        self.b_xyz = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_xyz,
                select='file',
                ftypes=[
                    ('Comma-Separated Values (.csv)', '*.csv'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_xyz.grid(
            sticky=W,
            row=0,
            column=2,
            pady=pad,
        )

        self.l_show = ttk.Label(
            root,
            text='Plot elevation profile:',
        )
        self.l_show.grid(
            stick=E,
            row=1,
            column=0,
            pady=pad,
        )
        self.e_show = ttk.Button(
            root,
            text='Plot!',
            command=lambda: show_plot(self.e_xyz.get()),
        )
        self.e_show.grid(
            stick=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.l_breaks = ttk.Label(
            root,
            text='Breakpoints (comma separated, no spaces)',
        )
        self.l_breaks.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )

        self.e_breaks = ttk.Entry(root)
        self.e_breaks.grid(
            stick=E,
            row=2,
            column=1,
            pady=pad,
        )
        self.e_breaks.insert(END, '')

        self.l_show = ttk.Label(
            root,
            text='Plot fit:',
        )
        self.l_show.grid(
            stick=E,
            row=3,
            column=0,
            pady=pad,
        )
        self.e_show = ttk.Button(
            root,
            text='Plot!',
            command=lambda: show_fit_plots(
                self.e_xyz.get(),
                self.e_breaks.get(),
            ),
        )
        self.e_show.grid(
            stick=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.l_dem = ttk.Label(
            root,
            text='DEM location:',
        )
        self.l_dem.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_dem2 = ttk.Entry(root)
        self.e_dem2.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )
        self.e_dem2.insert(END, '')
        self.e_dem2.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )

        self.b_dem2 = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_dem2,
                select='file',
                ftypes=[
                    ('TIFF, .tif', '*.tif'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_dem2.grid(
            sticky=W,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_clip = ttk.Label(
            root,
            text='DEM clip AOI (optional):',
        )
        self.l_clip.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_clip = ttk.Entry(root)
        self.e_clip.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
        )
        self.e_clip.insert(END, '')
        self.e_clip.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
            padx=5,
        )
        self.b_clip = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_clip,
                select='file',
                ftypes=[
                    ('Shapefile (.shp), .shp', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_clip.grid(
            sticky=W,
            row=5,
            column=2,
            pady=pad,
        )

        self.e_detrend = ttk.Button(
            root,
            text='Detrend DEM!',
            command=lambda: detrend(
                xyz_csv=self.e_xyz.get(),
                in_dem=self.e_dem2.get(),
                aoi_shp=self.e_clip.get(),
            ),
        )
        self.e_detrend.grid(
            sticky=E,
            row=6,
            column=0,
            pady=15,
        )
        root.grid_rowconfigure(
            6,
            minsize=50,
        )

        # Flow-stage modeling
        ######################################################################

        root = self.tabs[
            'Flow-stage modeling',
        ]

        self.top_label = ttk.Label(
            root,
            text='Run',
        )

        class WettedController:
            self.maxs = [0]
            self.runs = 0
            self.highest = max(self.maxs)
            self.switch = False

            def run(
                    self,
                    detrended_dem,
                    max_stage,
            ) -> None:
                self.runs += 1
                if max_stage > self.highest:
                    self.switch = True

                out_dir = os.path.dirname(detrended_dem)
                wet_dir = out_dir + '\\wetted_polygons'

                # ff you run for the first time it starts from 0 to max stage, generating shapefiles and plotting
                # or the updated max_stage is higher, shapefiles are made for the missing stages, and plots are updated
                if self.runs == 0 or self.switch:
                    print('Making wetted area polygons...')
                    prep_small_inc(
                        detrended_dem=detrended_dem,
                        max_stage=max_stage,
                    )
                    print('Done')

                    print('Creating flow stage analysis plots...')
                    imgs = pdf_cdf_plotting(
                        in_dir=wet_dir,
                        out_folder=out_dir,
                        max_stage=max_stage,
                    )
                    print('Done')
                    print('Plots @ %s' % os.path.dirname(imgs[0]))

                    for img in imgs:
                        name = os.path.basename(img)[:-4]
                        open_popup(name, img)
                        self.switch = False

                # ff the max stage is less than on equal to the previous max stage, plots are updated
                else:
                    print('Creating flow stage analysis plots...')
                    imgs = pdf_cdf_plotting(
                        in_dir=wet_dir,
                        out_folder=out_dir,
                        max_stage=max_stage,
                    )
                    print('Done')
                    print('Plots @ %s' % os.path.dirname(imgs[0]))

                    for img in imgs:
                        open_popup('Flow stage vs wetted area', img)

                # update list of used max stage values
                self.maxs.append(max_stage)

        # build GUI for flow-stage analysis
        self.l_detrended = ttk.Label(
            root,
            text='Detrended DEM:',
        )
        self.l_detrended.grid(
            sticky=E,
            row=0,
            column=0,
            pady=pad,
        )

        self.e_detrended = ttk.Entry(root)
        self.e_detrended.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_detrended.insert(END, '')

        self.b_detrended = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_detrended,
                select='file',
                ftypes=[
                    ('TIFF, .tif', '*.tif'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_detrended.grid(
            sticky=W,
            row=0,
            column=2,
            pady=pad,
        )

        self.l_max = ttk.Label(
            root,
            text='Max stage height:',
        )
        self.l_max.grid(
            sticky=E,
            row=1,
            column=0,
            pady=pad,
        )

        self.e_max = ttk.Entry(root)
        self.e_max.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_max.insert(END, 0)

        self.n_max = ttk.Label(
            root,
            text='Integer only, in DEM units',
        )
        self.n_max.grid(
            sticky=W,
            row=1,
            column=2,
            pady=pad,
        )

        self.e_flows = ttk.Button(
            root,
            text='Flow-stage analysis!',
            command=lambda: WettedController.run(
                self,
                detrended_dem=self.e_detrended.get(),
                max_stage=int(self.e_max.get()),
            ),
        )
        self.e_flows.grid(
            sticky=E,
            row=2,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(
            2,
            minsize=50,
        )

        self.note1 = ttk.Label(
            root,
            text='Choose key flow stages from plots and wetted area polygons',
        )
        self.note1.grid(
            sticky=W,
            row=3,
            columnspan=3,
            pady=pad,
        )

        self.l_zs = ttk.Label(
            root,
            text='Key stage heights:',
        )
        self.l_zs.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_zs = ttk.Entry(root)
        self.e_zs.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_zs.insert(END, '')

        self.n_zs = ttk.Label(
            root,
            text='Float only, comma separated, DEM units (ex: 0.6,1.7,5.8)',
        )
        self.n_zs.grid(
            sticky=E,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_dcenter = ttk.Label(
            root,
            text='Generate draft center-lines:',
        )
        self.l_dcenter.grid(
            stick=E,
            row=5,
            column=0,
            pady=15,
        )

        self.e_dcenter = ttk.Button(
            root,
            text='Run',
            command=lambda: stage_centerlines(
                dem=self.e_detrended.get(),
                zs=self.e_zs.get(),
                drafting=True,
            ),
        )
        self.e_dcenter.grid(
            sticky=E,
            row=5,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(
            2,
            minsize=50,
        )

        self.note2 = ttk.Label(
            root,
            text='Edit drafts center-lines with ArcGIS, then run below',
        )
        self.note2.grid(
            sticky=W,
            row=6,
            columnspan=3,
            pady=pad,
        )

        self.l_center = ttk.Label(
            root,
            text='Generate final center-lines:',
        )
        self.l_center.grid(
            stick=E,
            row=7,
            column=0,
            pady=15,
        )

        self.e_center = ttk.Button(
            root,
            text='Run',
            command=lambda: stage_centerlines(
                dem=self.e_detrended.get(),
                zs=self.e_zs.get(),
                drafting=False,
            ),
        )
        self.e_center.grid(
            sticky=E,
            row=7,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(
            2,
            minsize=50,
        )

        # Generate GCS series .csv files w/ landform classifications
        ######################################################################

        root = self.tabs[
            'GCS analysis',
        ]

        def gcs_analysis(
                detrended_dem,
                zs,
                xs_lengths,
                xs_spacing,
                analysis,
                clip_poly='',
                stage_plots=False,
                nest_plots=False,
        ) -> None:
            """DUMMY FUNCTION FOR FORMATTING"""
            if not analysis:
                print('Extract GCS series...')
                extract_gcs(
                    detrended_dem,
                    zs,
                    xs_lengths,
                    xs_spacing,
                    clip_poly=clip_poly,
                )
                print('Done')

            # TODO: finish this
            elif stage_plots and not nest_plots:
                print('Stage plots')
            elif stage_plots and nest_plots:
                print('Both plots')
            elif not stage_plots and nest_plots:
                print('Nest plots')

            print(stage_plots)
            print(nest_plots)
            print('In the gcs function')

        self.l_detrended2 = ttk.Label(
            root,
            text='Detrended DEM:',
        )
        self.l_detrended2.grid(
            sticky=E,
            row=0,
            column=0,
            pady=pad,
        )

        self.e_detrended2 = ttk.Entry(root)
        self.e_detrended2.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_detrended2.insert(END, '')

        self.b_detrended2 = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_detrended2,
                select='file',
                ftypes=[('TIFF, .tif', '*.tif'),
                        ('All files', '*'),
                        ],
            ),
        )
        self.b_detrended2.grid(
            sticky=W,
            row=0,
            column=2,
            pady=pad,
        )

        self.l_zs2 = ttk.Label(
            root,
            text='Key stage heights:',
        )
        self.l_zs2.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )
        self.e_zs2 = ttk.Entry(root)
        self.e_zs2.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_zs2.insert(END, '')

        self.n_zs2 = ttk.Label(
            root,
            text='Float only, comma separated, DEM units (ex: 0.6,1.7,5.8)',
        )
        self.n_zs2.grid(
            sticky=W,
            row=2,
            column=2,
            pady=pad,
        )

        self.l_length = ttk.Label(
            root,
            text='Cross-section lengths:',
        )
        self.l_length.grid(
            sticky=E,
            row=3,
            column=0,
            pady=pad,
        )

        self.e_length = ttk.Entry(root)
        self.e_length.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_length.insert(END, '')

        self.n_length = ttk.Label(
            root,
            text='List of integers corresponding to key stage heights (ex: 400,600,1000)',
        )
        self.n_length.grid(
            sticky=W,
            row=3,
            column=2,
            pady=pad,
        )

        self.l_space = ttk.Label(
            root,
            text='Cross-section spacing:',
        )
        self.l_space.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_space = ttk.Entry(root)
        self.e_space.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_space.insert(END, '')

        self.n_space = ttk.Label(
            root,
            text='Integer, in same units as the DEM. Should not be less than the DEM resolution!',
        )
        self.n_space.grid(
            sticky=W,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_clip = ttk.Label(
            root,
            text='Clip polygon (optional):',
        )
        self.l_clip.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_clip = ttk.Entry(root)
        self.e_clip.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_clip.insert(END, '')

        self.b_clip = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_clip,
                select='file',
                ftypes=[
                    ('Shapefiles', '*.shp'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_clip.grid(
            sticky=W,
            row=5,
            column=2,
            pady=pad,
        )

        self.n_clip = ttk.Label(
            root,
            text='Allows AOI to be updated, if unchanged leave empty!',
        )
        self.n_clip.grid(
            sticky=W,
            row=5,
            column=3,
            pady=pad,
        )

        self.l_gcs = ttk.Label(
            root,
            text='Extract GCS series:',
        )
        self.l_gcs.grid(
            stick=E,
            row=6,
            column=0,
            pady=15,
        )

        self.e_gcs = ttk.Button(
            root,
            text='Run',
            command=lambda: gcs_analysis(
                detrended_dem=self.e_detrended2.get(),
                zs=self.e_zs2.get(),
                xs_lengths=self.e_length.get(),
                xs_spacing=self.e_space.get(),
                clip_poly=self.e_clip.get(),
                analysis=False,
            ),
        )
        self.e_gcs.grid(
            sticky=E,
            row=6,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(
            6,
            minsize=50,
        )

        self.note2 = ttk.Label(
            root,
            text='Verify that cross-section lengths are sufficient before continuing! \
                Re-run above if necessary. Must have extracted GCS series first.',
        )
        self.note2.grid(
            sticky=W,
            row=7,
            columnspan=3,
            pady=pad,
        )

        self.l_plots = ttk.Label(
            root,
            text='Run GCS stage analysis?:',
        )
        self.l_plots.grid(
            sticky=E,
            row=8,
            column=0,
        )

        self.plots = BooleanVar()
        self.plots.set(False)

        self.r_plots_y = ttk.Radiobutton(
            root,
            text='Yes',
            variable=self.plots,
            value=True,
        )
        self.r_plots_y.grid(
            sticky=W,
            row=8,
            column=1,
        )

        self.r_plots_n = ttk.Radiobutton(
            root,
            text='No',
            variable=self.plots,
            value=False,
        )
        self.r_plots_n.grid(
            sticky=W,
            row=8,
            column=2,
            pady=pad,
        )
        root.grid_rowconfigure(
            8,
            minsize=30,
        )

        self.l_plots2 = ttk.Label(
            root,
            text='Run GCS nesting analysis?:',
        )
        self.l_plots2.grid(
            sticky=E,
            row=9,
            column=0,
        )

        self.plots2 = BooleanVar()
        self.plots2.set(False)

        self.r_plots_y2 = ttk.Radiobutton(
            root,
            text='Yes',
            variable=self.plots2,
            value=True,
        )
        self.r_plots_y2.grid(
            sticky=W,
            row=9,
            column=1,
        )

        self.r_plots_n2 = ttk.Radiobutton(
            root,
            text='No',
            variable=self.plots2,
            value=False,
        )
        self.r_plots_n2.grid(
            sticky=W,
            row=9,
            column=2,
            pady=pad,
        )
        root.grid_rowconfigure(
            9,
            minsize=30,
        )

        self.l_gcs = ttk.Label(
            root,
            text='GCS analysis:',
        )
        self.l_gcs.grid(
            stick=E,
            row=10,
            column=0,
            pady=15,
        )

        self.e_gcs = ttk.Button(
            root,
            text='Run',
            command=lambda: gcs_analysis(
                detrended_dem=self.e_detrended2.get(),
                zs=self.e_zs2.get(),
                xs_lengths=self.e_length.get(),
                xs_spacing=self.e_space.get(),
                analysis=False,
                clip_poly=self.e_clip.get(),
                stage_plots=self.plots.get(),
                nest_plots=self.plots2.get(),
            ),
        )
        self.e_gcs.grid(
            sticky=E,
            row=10,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(
            10,
            minsize=50,
        )

        # Generate river builder inputs from harmonic decomposition
        ######################################################################
        root = self.tabs[
            'River Builder prep',
        ]

        # TODO: make this work

        def river_builder_harmonics(
                in_csv,
                index_field,
                units,
                field_names,
                r_2,
                n,
                methods,
        ) -> None:
            """DUMMY FUNCTION FOR FORMATTING"""
            print('In the RB function')

        self.l_csv = ttk.Label(
            root,
            text='In csv:',
        )
        self.l_csv.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
        )

        self.e_csv = ttk.Entry(root)
        self.e_csv.insert(END, '')
        self.e_csv.grid(
            row=0,
            column=2,
            pady=pad,
        )

        self.b_csv = ttk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.b_csv,
                select='file',
                ftypes=[
                    ('Comma-delimited text', '*.csv'),
                    ('All files', '*'),
                ],
            ),
        )
        self.b_csv.grid(
            sticky=W,
            row=0,
            column=3,
            pady=pad,
        )

        self.l_field = ttk.Label(
            root,
            text='Index field:',
        )
        self.l_field.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.e_field = ttk.Entry(root)
        self.e_field.insert(END, '')
        self.e_field.grid(
            row=1,
            column=2,
            pady=pad,
        )

        self.l_units = ttk.Label(
            root,
            text='   Units:',
        )
        self.l_units.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
        )

        self.e_units = StringVar()

        self.r_meters = ttk.Radiobutton(
            root,
            text='Meters',
            variable=self.e_units,
            value='m',
        )
        self.r_meters.grid(
            row=2,
            column=2,
            pady=pad,
        )
        self.r_feet = ttk.Radiobutton(
            root,
            text='US Feet',
            variable=self.e_units,
            value='ft',
        )
        self.r_feet.grid(
            row=2,
            column=3,
            pady=pad,
        )

        self.l_labels = ttk.Label(
            root,
            text='Add signal labels (optional, comma separated!)',
        )
        self.l_labels.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.e_labels = ttk.Entry(root)
        self.e_labels.insert(END, '')
        self.e_labels.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.l_r2 = ttk.Label(
            root,
            text='R^2 threshold:',
        )
        self.l_r2.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )

        self.e_r2 = ttk.Entry(root)
        self.e_r2.insert(END, 0.90)
        self.e_r2.grid(
            row=4,
            column=2,
            pady=pad,
        )

        self.l_harms = ttk.Label(
            root,
            text='N harmonics override (optional, leave at 0):',
        )
        self.l_harms.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
        )

        self.e_harms = ttk.Entry(root)
        self.e_harms.insert(END, 0)
        self.e_harms.grid(
            row=5,
            column=2,
            pady=pad,
        )

        self.l_meth = ttk.Label(
            root,
            text='Select interpolation method:',
        )
        self.l_meth.grid(
            sticky=E,
            row=6,
            column=2,
            pady=pad,
        )

        # linear and natural neighbors refer to TIN based methods, be sure to document
        methods2 = [
            'by_fft',
            'by_power',
            'by_power_binned',
        ]

        self.meth = StringVar()

        self.e_meth = ttk.OptionMenu(
            root,
            self.meth,
            *methods2,
        )
        self.e_meth.grid(
            sticky=W,
            row=6,
            column=3,
            pady=pad,
        )

        b = Button(
            root,
            text='   Run    ',
            command=lambda: river_builder_harmonics(
                in_csv=str.replace(self.e_csv.get(), "\\", "\\\\"),
                index_field=self.e_field.get(),
                units=self.e_units.get(),
                field_names=string_to_list(self.e_labels.get()),
                r_2=float(self.e_r2.get()),
                n=int(self.e_harms.get()),
                methods=self.meth.get(),
            ),
        )
        b.grid(
            sticky=W,
            row=7,
            column=2,
        )
        root.grid_rowconfigure(
            4,
            minsize=80,
        )


if __name__ == '__main__':
    GCSGraphicUserInterface().mainloop()
