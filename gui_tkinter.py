import os
import tkinter as tk
from tkinter import filedialog
from typing import List, Tuple
from PIL import Image, ImageTk

from gcs_analysis_tools.utils import init_logger

from gcs_analysis_tools.lidar_to_dem.gui_funcs import lidar_prep, dem_generation
#from gcs_analysis_tools.detrend_dem.thalweg_creation import *
from gcs_analysis_tools.detrend_dem.gui_funcs import make_xyz_plot, make_fit_plots, detrend
from gcs_analysis_tools.flow_stages.gui_funcs import WetterController, model_each_flow_stage, stage_centerlines
from gcs_analysis_tools.gcs_analysis.gui_funcs import run_gcs_analyses
from gcs_analysis_tools.river_builder_prep.gui_funcs import export_to_river_builder 


class GCSGraphicUserInterface(tk.Frame):

    # Initialize the class, this allows the GUI to run when the code is ran
    def __init__(
            self,
            master=None,
    ) -> None:

        # Initialize the tk frame that will hold all tabs holding each processing step
        tk.Frame.__init__(
            self,
            master,
        )
        self.pack()

        # initialize attribute files CHANGE
        self.dem = tk.StringVar()
        self.det_dem = tk.StringVar()
        self.centerline = tk.StringVar()
        self.station_lines = tk.StringVar()

        # set header
        self.master.title('Geomorphic Covariance Structure (GCS) analysis GUI')
        self.master.iconbitmap('imgs/win_icon.ico')

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
        self.style = tk.Style()

        # Adding the Breeze tk theme https://github.com/MaxPerl/ttk-Breeze
        breeze_dir = os.getcwd() + '\\tk-Breeze-master'

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
        self.tab_container = tk.Notebook(master)

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
            tab = tk.Frame(self.tab_container)
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

            top = tk.Toplevel(root)
            top.geometry()
            top.title(title)

            self.ph = ImageTk.PhotoImage(self.im, master=top)
            self.label = tk.Label(top, image=self.ph)
            self.label.image = self.ph
            self.label.grid(row=1, column=1, columnspan=3)

            self.label2 = tk.Label(top, text='Image saved @ %s' % image)
            self.label2.grid(row=2, column=1)

        # LiDAR prep (filling tabs w/ widgets)
        ######################################################################

        root = self.tabs['LiDAR Data prep']

        self.l_lasbin1 = tk.Label(
            root,
            text='LAStools /bin/ directory:',
        )
        self.l_lasbin1.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
        )

        self.e_lasbin1 = tk.Entry(root)
        self.e_lasbin1.insert(END, str(os.getcwd() + '\\LAStools\\bin'))
        self.e_lasbin1.grid(
            row=0,
            column=2,
            pady=pad,
        )

        self.b_lasbin1 = tk.Button(
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

        self.l_lidardir1 = tk.Label(
            root,
            text='LiDAR data directory:',
        )
        self.l_lidardir1.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.e_lidardir1 = tk.Entry(root)
        self.e_lidardir1.insert(END, '')
        self.e_lidardir1.grid(
            row=1,
            column=2,
            pady=pad,
        )
        self.b_lidardir1 = tk.Button(
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

        self.l_in_spatialref = tk.Label(
            root,
            text='LiDAR spatial reference (.shp):',
        )
        self.l_in_spatialref.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
        )

        self.e_in_spatialref = tk.Entry(root)
        self.e_in_spatialref.insert(END, '')
        self.e_in_spatialref.grid(
            row=2,
            column=2,
            pady=pad,
        )

        self.b_in_spatialref = tk.Button(
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

        self.l_naip = tk.Label(
            root,
            text='NAIP imagery folder:',
        )
        self.l_naip.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.e_naip = tk.Entry(root)
        self.e_naip.insert(END, '')
        self.e_naip.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.b_naip = tk.Button(
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

        self.l_ndvi = tk.Label(
            root,
            text='NDVI vegetation threshold:',
        )
        self.l_ndvi.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )
        self.e_ndvi = tk.Entry(root)
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

        self.l_aoi = tk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_aoi.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
        )

        self.e_aoi = tk.Entry(root)
        self.e_aoi.insert(END, '')
        self.e_aoi.grid(
            row=5,
            column=2,
            pady=pad,
        )

        self.b_aoi = tk.Button(
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

        self.spacer1 = tk.Label(root, text='')
        self.spacer1.grid(
            sticky=W,
            row=6,
            column=1,
            pady=pad,
        )
        self.prep_run = tk.Button(
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
        root.grid_rowconfigure(17, minsize=80)

        self.instruct = tk.Label(
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
        self.l_lasbin = tk.Label(
            root,
            text='LAStools /bin/ directory:',
        )
        self.l_lasbin.grid(
            sticky=E,
            row=0,
            column=1,
        )

        self.e_lasbin = tk.Entry(root)
        self.e_lasbin.insert(END, os.getcwd() + '\\LAStools\\bin')
        self.e_lasbin.grid(
            row=0,
            column=2,
        )
        self.b_lasbin = tk.Button(
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

        self.l_lidardir = tk.Label(
            root,
            text='LiDAR data directory:',
        )
        self.l_lidardir.grid(
            sticky=E,
            row=1,
            column=1,
        )

        self.e_lidardir = tk.Entry(root)
        self.e_lidardir.insert(END, '')
        self.e_lidardir.grid(
            row=1,
            column=2,
        )

        self.b_lidardir = tk.Button(
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

        self.shp_var = tk.StringVar()

        self.l_ground_shp = tk.Label(
            root,
            text='Ground polygon (.shp):',
        )
        self.l_ground_shp.grid(
            sticky=E,
            row=2,
            column=1,
        )

        self.e_ground_shp = tk.Entry(
            root,
            textvariable=self.shp_var,
        )
        self.e_ground_shp.grid(
            row=2,
            column=2,
        )

        self.b_ground_shp = tk.Button(
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

        self.l_out_spatialref = tk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_out_spatialref.grid(
            sticky=E,
            row=3,
            column=1,
        )

        self.e_out_spatialref = tk.Entry(root)
        self.e_out_spatialref.insert(END, '')
        self.e_out_spatialref.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.b_out_spatialref = tk.Button(
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

        self.l_coarse_class = tk.Label(
            root,
            text='standard/coarse classification parameters:',
        )
        self.l_coarse_class.grid(
            row=5,
            column=0,
            columnspan=2,
        )

        self.l_c_step = tk.Label(
            root,
            text='step size:',
        )
        self.l_c_step.grid(
            sticky=E,
            row=6,
        )

        self.e_c_step = tk.Entry(root)
        self.e_c_step.grid(
            row=6,
            column=1,
        )

        self.l_c_bulge = tk.Label(
            root,
            text='bulge:',
        )
        self.l_c_bulge.grid(
            sticky=E,
            row=7,
        )

        self.e_c_bulge = tk.Entry(root)
        self.e_c_bulge.grid(
            row=7,
            column=1,
        )

        self.l_c_spike = tk.Label(
            root,
            text='spike:',
        )
        self.l_c_spike.grid(
            sticky=E,
            row=8,
        )
        self.e_c_spike = tk.Entry(root)
        self.e_c_spike.grid(
            row=8,
            column=1,
        )

        self.l_c_dspike = tk.Label(
            root,
            text='down spike:',
        )
        self.l_c_dspike.grid(
            sticky=E,
            row=9,
        )

        self.e_c_dspike = tk.Entry(root)
        self.e_c_dspike.grid(
            row=9,
            column=1,
        )

        self.l_c_offset = tk.Label(
            root,
            text='offset:',
        )
        self.l_c_offset.grid(
            sticky=E,
            row=10,
        )
        self.e_c_offset = tk.Entry(root)
        self.e_c_offset.grid(
            row=10,
            column=1,
        )

        self.l_fine_class = tk.Label(
            root,
            text='fine classification parameters (in ground area):',
        )
        self.l_fine_class.grid(
            row=5,
            column=2,
            columnspan=2,
        )

        self.l_f_step = tk.Label(
            root,
            text='step size:',
        )
        self.l_f_step.grid(
            sticky=E,
            row=6,
            column=2,
        )

        self.e_f_step = tk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_step.grid(
            row=6,
            column=3,
        )

        self.l_f_bulge = tk.Label(
            root,
            text='bulge:',
        )
        self.l_f_bulge.grid(
            sticky=E,
            row=7,
            column=2,
        )

        self.e_f_bulge = tk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_bulge.grid(
            row=7,
            column=3,
        )

        self.l_f_spike = tk.Label(
            root,
            text='spike:',
        )
        self.l_f_spike.grid(
            sticky=E,
            row=8,
            column=2,
        )

        self.e_f_spike = tk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_spike.grid(
            row=8,
            column=3,
        )

        self.l_f_dspike = tk.Label(
            root,
            text='down spike:',
        )
        self.l_f_dspike.grid(
            sticky=E,
            row=9,
            column=2,
        )

        self.e_f_dspike = tk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_dspike.grid(
            row=9,
            column=3,
        )

        self.l_f_offset = tk.Label(
            root,
            text='offset:',
        )
        self.l_f_offset.grid(
            sticky=E,
            row=10,
            column=2,
        )

        self.e_f_offset = tk.Entry(
            root,
            state=DISABLED,
        )
        self.e_f_offset.grid(
            row=10,
            column=3,
        )

        # specify units
        self.l_lidar_units = tk.Label(
            root,
            text='Units:',
        )
        self.l_lidar_units.grid(
            sticky=W,
            row=11,
            column=2,
        )
        root.grid_rowconfigure(11, minsize=30)

        self.lidar_units = tk.StringVar()

        self.r_lidar_meters = tk.Radiobutton(
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
        self.r_lidar_feet = tk.Radiobutton(
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
        self.l_lidar_cores = tk.Label(
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

        self.r1_lidar = tk.Radiobutton(
            root,
            text='1',
            variable=self.core_num, value=1,
        )
        self.r1_lidar.grid(
            sticky=E,
            row=14,
            column=1,
        )

        self.r2_lidar = tk.Radiobutton(
            root,
            text='2',
            variable=self.core_num,
            value=2,
        )
        self.r2_lidar.grid(
            row=14,
            column=2,
        )

        self.r4_lidar = tk.Radiobutton(
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

        self.r8_lidar = tk.Radiobutton(
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

        self.r16_lidar = tk.Radiobutton(
            root,
            text='16',
            variable=self.core_num,
            value=16,
        )
        self.r16_lidar.grid(
            row=15,
            column=2,
        )

        self.r32_lidar = tk.Radiobutton(
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

        self.l_keep_orig_lidar = tk.Label(
            root,
            text='Keep original ground/veg points: ',
        )
        self.l_keep_orig_lidar.grid(
            sticky=E,
            row=16,
            column=1,
        )

        self.keep_orig_lidar = tk.BooleanVar()

        self.c_keep_orig_lidar = tk.Checkbutton(
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

        self.l_dem_res = tk.Label(
            root,
            text='DEM resolution (meters):',
        )
        self.l_dem_res.grid(
            sticky=E,
            row=17,
            column=2,
        )

        self.e_dem_res = tk.Entry(root)
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

        self.l_dem_meth = tk.Label(
            root,
            text='Select interpolation method:',
        )
        self.l_dem_meth.grid(
            sticky=E,
            row=18,
            column=1,
            pady=pad,
        )

        self.e_dem_meth = tk.StringVar()

        self.option_menu1 = tk.OptionMenu(
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

        self.l_void_meth = tk.Label(
            root,
            text='Void fill method (for binning:',
        )
        self.l_void_meth.grid(
            sticky=E,
            row=19,
            column=1,
            pady=pad,
        )

        self.e_void_meth = tk.StringVar()

        self.option_menu2 = tk.OptionMenu(
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

        self.l_tri_meth = tk.Label(
            root,
            text='Triangulation method:',
        )
        self.l_tri_meth.grid(
            sticky=E,
            row=20,
            column=1,
            pady=pad,
        )
        self.e_tri_meth = tk.StringVar()

        self.option_menu3 = tk.OptionMenu(
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

        # make 'Run' tk.Button in GUI to call the process_lidar() function
        self.b_lidar_run = tk.Button(
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
        root.grid_rowconfigure(20, minsize=40)

        # Generate thalweg centerline and extract elevation profile
        ######################################################################
        root = self.tabs['Thalweg centerline']

        self.remind = tk.Label(
            root,
            text='Create upstream flow polygon in ArcMap/Pro',
        )
        self.remind.grid(
            sticky=E,
            row=0,
            column=0,
        )

        self.l_flow_poly = tk.Label(
            root,
            text='Upstream flow polygon (.shp):',
        )
        self.l_flow_poly.grid(
            sticky=E,
            row=1,
            column=0,
            pady=pad,
        )

        self.e_flow_poly = tk.Entry(root)
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

        self.b_flow_poly = tk.Button(
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

        self.l_extent = tk.Label(
            root,
            text='AOI shapefile (.shp):',
        )
        self.l_extent.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )

        self.e_extent = tk.Entry(root)
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

        self.b_extent = tk.Button(
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

        self.l_filt = tk.Label(
            root,
            text='Filter passes (15x default):',
        )
        self.l_filt.grid(
            sticky=E,
            row=3,
            column=0,
            pady=pad,
        )

        self.e_filt = tk.Entry(root)
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

        self.l_smooth = tk.Label(
            root,
            text='Smoothing distance (meters):',
        )
        self.l_smooth.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_smooth = tk.Entry(root)
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

        self.l_dem = tk.Label(
            root,
            text='DEM (.tif):',
        )
        self.l_dem.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_dem = tk.Entry(root)
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

        self.b_dem = tk.Button(
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
        self.b_detrend_prep1 = tk.Button(
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

        self.l_step = tk.Label(
            root,
            text='Verify centerline quality (edit if necessary), then run below...',
        )
        self.l_step.grid(
            sticky=E,
            row=7,
            column=0,
        )

        # Create run button to generate a thalweg elevation table from 1m spaced station points
        self.b_detrend_prep2 = tk.Button(
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
        root.grid_rowconfigure(6, minsize=50)

        # DEM detrending
        ######################################################################
        root = self.tabs['Detrend DEM']
        self.l_xyz = tk.Label(
            root,
            text='Thalweg profile csv:',
        )
        self.l_xyz.grid(
            sticky=E,
            row=0,
            column=0,
        )

        self.e_xyz = tk.Entry(root)
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

        self.b_xyz = tk.Button(
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

        self.l_show = tk.Label(
            root,
            text='Plot elevation profile:',
        )
        self.l_show.grid(
            stick=E,
            row=1,
            column=0,
            pady=pad,
        )
        self.e_show = tk.Button(
            root,
            text='Plot!',
            command=lambda: open_popup('Thalweg elevation profile', make_xyz_plot(self.e_xyz.get())),
        )
        self.e_show.grid(
            stick=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.l_breaks = tk.Label(
            root,
            text='Breakpoints (comma separated, no spaces)',
        )
        self.l_breaks.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )

        self.e_breaks = tk.Entry(root)
        self.e_breaks.grid(
            stick=E,
            row=2,
            column=1,
            pady=pad,
        )
        self.e_breaks.insert(END, '')

        self.l_show = tk.Label(
            root,
            text='Plot fit:',
        )
        self.l_show.grid(
            stick=E,
            row=3,
            column=0,
            pady=pad,
        )

        def show_fit_plots(fit_plot, res_plot, breakpoint_list):
            open_popup(
                'Linear fit w/ breakpoints: %s' % breakpoint_list,
                fit_plot,
            )
            open_popup(
                'Residual plot w/ breakpoints: %s' % breakpoint_list,
                res_plot,
            )
        self.e_show = tk.Button(
            root,
            text='Plot!',
            command=lambda: show_fit_plots(make_fit_plots(
                self.e_xyz.get(),
                self.e_breaks.get(),
            )),
        )
        self.e_show.grid(
            stick=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.l_dem = tk.Label(
            root,
            text='DEM location:',
        )
        self.l_dem.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_dem2 = tk.Entry(root)
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

        self.b_dem2 = tk.Button(
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

        self.l_clip = tk.Label(
            root,
            text='DEM clip AOI (optional):',
        )
        self.l_clip.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_clip = tk.Entry(root)
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
        self.b_clip = tk.Button(
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

        self.e_detrend = tk.Button(
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
        root.grid_rowconfigure(6, minsize=50)

        # Flow-stage modeling
        ######################################################################

        root = self.tabs['Flow-stage modeling']

        self.top_label = tk.Label(
            root,
            text='Run',
        )
        # build GUI for flow-stage analysis
        self.l_detrended = tk.Label(
            root,
            text='Detrended DEM:',
        )
        self.l_detrended.grid(
            sticky=E,
            row=0,
            column=0,
            pady=pad,
        )

        self.e_detrended = tk.Entry(root)
        self.e_detrended.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_detrended.insert(END, '')

        self.b_detrended = tk.Button(
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

        self.l_max = tk.Label(
            root,
            text='Max stage height:',
        )
        self.l_max.grid(
            sticky=E,
            row=1,
            column=0,
            pady=pad,
        )

        self.e_max = tk.Entry(root)
        self.e_max.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_max.insert(END, 0)

        self.n_max = tk.Label(
            root,
            text='Integer only, in DEM units',
        )
        self.n_max.grid(
            sticky=W,
            row=1,
            column=2,
            pady=pad,
        )
        controller = WetterController()
        def show_flow_stage_plots(named_imgs: List[Tuple[str, str]]):
            # TODO: combine this with the other simular func
            for name, img in named_imgs:
                open_popup(name, img)

        self.e_flows = tk.Button(
            root,
            text='Flow-stage analysis!',
            command=lambda: show_flow_stage_plots(model_each_flow_stage(
                detrended_dem=self.e_detrended.get(),
                max_stage=int(self.e_max.get()),
                controller=controller,
            )),
        )
        self.e_flows.grid(
            sticky=E,
            row=2,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(2, minsize=50)

        self.note1 = tk.Label(
            root,
            text='Choose key flow stages from plots and wetted area polygons',
        )
        self.note1.grid(
            sticky=W,
            row=3,
            columnspan=3,
            pady=pad,
        )

        self.l_zs = tk.Label(
            root,
            text='Key stage heights:',
        )
        self.l_zs.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_zs = tk.Entry(root)
        self.e_zs.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_zs.insert(END, '')

        self.n_zs = tk.Label(
            root,
            text='Float only, comma separated, DEM units (ex: 0.6,1.7,5.8)',
        )
        self.n_zs.grid(
            sticky=E,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_dcenter = tk.Label(
            root,
            text='Generate draft center-lines:',
        )
        self.l_dcenter.grid(
            stick=E,
            row=5,
            column=0,
            pady=15,
        )

        self.e_dcenter = tk.Button(
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
        root.grid_rowconfigure(2, minsize=50)

        self.note2 = tk.Label(
            root,
            text='Edit drafts center-lines with ArcGIS, then run below',
        )
        self.note2.grid(
            sticky=W,
            row=6,
            columnspan=3,
            pady=pad,
        )

        self.l_center = tk.Label(
            root,
            text='Generate final center-lines:',
        )
        self.l_center.grid(
            stick=E,
            row=7,
            column=0,
            pady=15,
        )

        self.e_center = tk.Button(
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
        root.grid_rowconfigure(2, minsize=50)

        # Generate GCS series .csv files w/ landform classifications
        ######################################################################

        root = self.tabs['GCS analysis']
        self.l_detrended2 = tk.Label(
            root,
            text='Detrended DEM:',
        )
        self.l_detrended2.grid(
            sticky=E,
            row=0,
            column=0,
            pady=pad,
        )

        self.e_detrended2 = tk.Entry(root)
        self.e_detrended2.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_detrended2.insert(END, '')

        self.b_detrended2 = tk.Button(
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

        self.l_zs2 = tk.Label(
            root,
            text='Key stage heights:',
        )
        self.l_zs2.grid(
            sticky=E,
            row=2,
            column=0,
            pady=pad,
        )
        self.e_zs2 = tk.Entry(root)
        self.e_zs2.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_zs2.insert(END, '')

        self.n_zs2 = tk.Label(
            root,
            text='Float only, comma separated, DEM units (ex: 0.6,1.7,5.8)',
        )
        self.n_zs2.grid(
            sticky=W,
            row=2,
            column=2,
            pady=pad,
        )

        self.l_length = tk.Label(
            root,
            text='Cross-section lengths:',
        )
        self.l_length.grid(
            sticky=E,
            row=3,
            column=0,
            pady=pad,
        )

        self.e_length = tk.Entry(root)
        self.e_length.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_length.insert(END, '')

        self.n_length = tk.Label(
            root,
            text='List of integers corresponding to key stage heights (ex: 400,600,1000)',
        )
        self.n_length.grid(
            sticky=W,
            row=3,
            column=2,
            pady=pad,
        )

        self.l_space = tk.Label(
            root,
            text='Cross-section spacing (integer):',
        )
        self.l_space.grid(
            sticky=E,
            row=4,
            column=0,
            pady=pad,
        )

        self.e_space = tk.Entry(root)
        self.e_space.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_space.insert(END, '')

        self.n_space = tk.Label(
            root,
            text='Integer, in same units as the DEM. Should not be less than the DEM resolution!',
        )
        self.n_space.grid(
            sticky=W,
            row=4,
            column=2,
            pady=pad,
        )

        self.l_clip = tk.Label(
            root,
            text='Clip polygon (optional):',
        )
        self.l_clip.grid(
            sticky=E,
            row=5,
            column=0,
            pady=pad,
        )

        self.e_clip = tk.Entry(root)
        self.e_clip.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
            padx=5,
        )
        self.e_clip.insert(END, '')

        self.b_clip = tk.Button(
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

        self.l_gcs = tk.Label(
            root,
            text='Extract GCS series:',
        )
        self.l_gcs.grid(
            stick=E,
            row=6,
            column=0,
            pady=15,
        )

        self.e_gcs = tk.Button(
            root,
            text='Run',
            command=lambda: run_gcs_analyses(
                detrended_dem=self.e_detrended2.get(),
                zs=self.e_zs2.get(),
                xs_lengths=self.e_length.get(),
                xs_spacing=self.e_space.get(),
                clip_poly=self.e_clip.get(),
            ),
        )
        self.e_gcs.grid(
            sticky=E,
            row=6,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(6, minsize=50)

        self.note2 = tk.Label(
            root,
            text=(
                'Verify that cross-section lengths are sufficient before continuing! '
                'Re-run above if necessary. Must have extracted GCS series first.'
            ),
        )
        self.note2.grid(
            sticky=W,
            row=7,
            columnspan=3,
            pady=pad,
        )

        self.l_plots = tk.Label(
            root,
            text='Run GCS stage analysis?:',
        )
        self.l_plots.grid(
            sticky=E,
            row=8,
            column=0,
        )

        self.plots = tk.BooleanVar()
        self.plots.set(False)

        self.r_plots_y = tk.Radiobutton(
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

        self.r_plots_n = tk.Radiobutton(
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
        root.grid_rowconfigure(8, minsize=30)

        self.l_plots2 = tk.Label(
            root,
            text='Run GCS nesting analysis?:',
        )
        self.l_plots2.grid(
            sticky=E,
            row=9,
            column=0,
        )

        self.plots2 = tk.BooleanVar()
        self.plots2.set(False)

        self.r_plots_y2 = tk.Radiobutton(
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

        self.r_plots_n2 = tk.Radiobutton(
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
        root.grid_rowconfigure(9, minsize=30)

        # choose where to put output analyses files
        self.l_analysis_dir = tk.Label(
            root,
            text='Override Analysis Output Directory (optional):',
        )
        self.l_analysis_dir.grid(
            stick=E,
            row=10,
            column=0,
            pady=15,
        )

        self.e_analysis_dir = tk.Entry(root)
        self.e_analysis_dir.insert(
            END,
            '',
        )
        self.e_analysis_dir.grid(
            row=10,
            column=1,
            pady=pad,
            padx=5,
        )

        self.b_analysis_dir = tk.Button(
            root,
            text='Browse',
            command=lambda: browse(
                root,
                self.e_analysis_dir,
                select='folder',
            ),
        )

        self.b_analysis_dir.grid(
            sticky=W,
            row=10,
            column=2,
            pady=pad,
        )

        # radio buttons to control which analyses to run
        self.l_gcs = tk.Label(
            root,
            text='GCS analysis:',
        )
        self.l_gcs.grid(
            stick=E,
            row=11,
            column=0,
            pady=15,
        )

        self.e_gcs = tk.Button(
            root,
            text='Run',
            command=lambda: run_gcs_analyses(
                detrended_dem=self.e_detrended2.get(),
                zs=self.e_zs2.get(),
                xs_lengths=self.e_length.get(),
                xs_spacing=self.e_space.get(),
                clip_poly=self.e_clip.get(),
                stage_plots=self.plots.get(),
                nest_plots=self.plots2.get(),
                analysis_dir=self.e_analysis_dir.get(),
            ),
        )
        self.e_gcs.grid(
            sticky=E,
            row=11,
            column=1,
            pady=15,
        )
        root.grid_rowconfigure(11, minsize=50)

        # Generate river builder inputs from harmonic decomposition
        ######################################################################
        root = self.tabs['River Builder prep']
        
        self.l_csv = tk.Label(
            root,
            text='In csv:',
        )
        self.l_csv.grid(
            sticky=E,
            row=0,
            column=1,
            pady=pad,
        )

        self.e_csv = tk.Entry(root)
        self.e_csv.insert(END, '')
        self.e_csv.grid(
            row=0,
            column=2,
            pady=pad,
        )

        self.b_csv = tk.Button(
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

        self.l_field = tk.Label(
            root,
            text='Index field:',
        )
        self.l_field.grid(
            sticky=E,
            row=1,
            column=1,
            pady=pad,
        )

        self.e_field = tk.Entry(root)
        self.e_field.insert(END, 'dist_down')
        self.e_field.grid(
            row=1,
            column=2,
            pady=pad,
        )

        self.l_units = tk.Label(
            root,
            text='   Units:',
        )
        self.l_units.grid(
            sticky=E,
            row=2,
            column=1,
            pady=pad,
        )

        self.e_units = tk.StringVar()

        self.r_meters = tk.Radiobutton(
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
        self.r_feet = tk.Radiobutton(
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

        self.l_labels = tk.Label(
            root,
            text='Add list of columns to export (comma separated, overrides W + Z))',
        )
        self.l_labels.grid(
            sticky=E,
            row=3,
            column=1,
            pady=pad,
        )

        self.e_labels = tk.Entry(root)
        self.e_labels.insert(END, '')
        self.e_labels.grid(
            row=3,
            column=2,
            pady=pad,
        )

        self.l_r2 = tk.Label(
            root,
            text='R^2 threshold:',
        )
        self.l_r2.grid(
            sticky=E,
            row=4,
            column=1,
            pady=pad,
        )

        self.e_r2 = tk.Entry(root)
        self.e_r2.insert(END, 0.90)
        self.e_r2.grid(
            row=4,
            column=2,
            pady=pad,
        )

        self.l_harms = tk.Label(
            root,
            text='N harmonics override (optional, leave at 0):',
        )
        self.l_harms.grid(
            sticky=E,
            row=5,
            column=1,
            pady=pad,
        )

        self.e_harms = tk.Entry(root)
        self.e_harms.insert(END, 0)
        self.e_harms.grid(
            row=5,
            column=2,
            pady=pad,
        )

        self.l_meth = tk.Label(
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

        self.meth = tk.StringVar()

        self.e_meth = tk.OptionMenu(
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

        b = tk.Button(
            root,
            text='   Run    ',
            command=lambda: export_to_river_builder(
                in_csv=str(self.e_csv.get()),
                index_field=self.e_field.get(),
                units=self.e_units.get(),
                r2_threshold=self.e_r2.get(),
                n_harmonics=self.e_harms.get(),
                methods=self.meth.get(),
                field_headers=self.e_labels.get(),
            ),
        )
        b.grid(
            sticky=W,
            row=7,
            column=2,
        )
        root.grid_rowconfigure(3, minsize=80)


if __name__ == '__main__':
    # initialize the logger
    init_logger(__file__)

    GCSGraphicUserInterface().mainloop()
