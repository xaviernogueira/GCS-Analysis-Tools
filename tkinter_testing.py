import tkinter as tk
from pydoc import browse
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
# from ttkthemes import ThemedTk
import warnings


# warnings.filterwarnings("ignore")


def browse(root, entry, select='file', ftypes=[('All files', '*')]):
    """GUI button command: opens browser window and adds selected file/folder to entry"""
    if select == 'file':
        filename = filedialog.askopenfilename(parent=root, title='Choose a file', filetypes=ftypes)
        if filename is not None:
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


# def spatial_ref(in_file):
# return arcpy.Describe(in_file).spatialReference

####### LETS MOVE THIS TO THE PROCESSING ITSELF, TOO COMPLICATED TO INSERT A SPATIAL REFERENCE OBJECT WITHIN A TKINTER BOX


class gcs_gui(tk.Frame):
    def __init__(self, master=None):  # Initialize the class, this allows the GUI to run when the code is ran
        tk.Frame.__init__(self, master)  # Initialize the tk frame that will hold all tabs holding each processing step
        self.pack()

        # initialize attribute files CHANGE
        self.dem = StringVar()
        self.det_dem = StringVar()
        self.centerline = StringVar()
        self.station_lines = StringVar()

        # set header
        self.master.title('Geomorphic Covariance Structure (GCS) analysis GUI')
        # self.master.iconbitmap('river.ico')  # UPDATE

        self.bg_color = 'LightSkyBlue1'  # Create color scheme and padding for the frame/self
        self.padding = 5

        ww = self.master.winfo_screenwidth() / 2  # controls window width
        wh = self.master.winfo_screenwidth() / 2  # controls window height
        wx = (self.master.winfo_screenwidth() - ww) / 2  # position relative to screen width and ww
        wy = (self.master.winfo_screenheight() - wh) / 2  # position relative to screen height and wh
        self.master.geometry("%dx%d+%d+%d" % (ww, wh, wx, wy))  # set window height and location

        # set widget styles
        self.style = ttk.Style()

        # Adding the Breeze ttk theme https://github.com/MaxPerl/ttk-Breeze
        self.tk.call('lappend', 'auto_path', r'C:\Users\xavie\Documents\ttk-Breeze-master')
        self.tk.call('source', r'C:\Users\xavie\Documents\ttk-Breeze-master\breeze.tcl')
        self.style.theme_use('Breeze')

        # initialize tab handler
        self.tab_container = ttk.Notebook(master)

        self.tab_names = ['LiDAR Data prep', 'DEM generation', 'Thalweg centerline',
                          'Detrend DEM', 'Flow-stage modeling', 'Extract GCS series and landforms',
                          'Flow-stage analysis', 'River Builder prep']

        self.tabs = {}
        for tab_name in self.tab_names:
            tab = ttk.Frame(self.tab_container)
            self.tab_container.add(tab, text=tab_name)
            self.tab_container.pack(expand=1, fill="both")
            self.tabs[tab_name] = tab

        pad = 5  # Denoted padding between tkinter widgets

        # LiDAR prep (filling tabs w/ widgets)
        ######################################################################

        root = self.tabs['LiDAR Data prep']

        def lidar_prep(lasbin, lidardir, spatial_shp, naip_folder, ndvi_thresh=0.4, aoi_shp=''):
            """This function is ran by the LiDAR Data prep tab.
             Outputs: A 1m (or other resolution) DEM modeling bare ground LiDAR returns """

            print('Running, functions go here...')

        self.l_lasbin = ttk.Label(root, text='LAStools /bin/ directory:')
        self.l_lasbin.grid(sticky=E, row=0, column=1, pady=pad)
        self.e_lasbin = ttk.Entry(root)
        self.e_lasbin.insert(END, '')
        self.e_lasbin.grid(row=0, column=2, pady=pad)
        self.b_lasbin = ttk.Button(root, text='Browse',
                                   command=lambda: browse(root, entry=self.e_lasbin, select='folder'))
        self.b_lasbin.grid(sticky=W, row=0, column=3, pady=pad)

        self.l_lidardir = ttk.Label(root, text='LiDAR data directory:')
        self.l_lidardir.grid(sticky=E, row=1, column=1, pady=pad)
        self.e_lidardir = ttk.Entry(root)
        self.e_lidardir.insert(END, '')
        self.e_lidardir.grid(row=1, column=2, pady=pad)
        self.b_lidardir = ttk.Button(root, text='Browse',
                                     command=lambda: browse(root, self.e_lidardir, select='folder'))
        self.b_lidardir.grid(sticky=W, row=1, column=3, pady=pad)

        self.l_in_spatialref = ttk.Label(root, text='LiDAR project spatial ref shapefile:')
        self.l_in_spatialref.grid(sticky=E, row=2, column=1, pady=pad)
        self.e_in_spatialref = ttk.Entry(root)
        self.e_in_spatialref.insert(END, '')
        self.e_in_spatialref.grid(row=2, column=2, pady=pad)
        self.b_in_spatialref = ttk.Button(root, text='Browse',
                                          command=lambda: browse(root, self.e_in_spatialref, select='file',
                                                                 ftypes=[('Shapefile', '*.shp'),
                                                                         ('All files', '*')]))
        self.b_in_spatialref.grid(sticky=W, row=2, column=3, pady=pad)

        self.l_naip = ttk.Label(root, text='NAIP imagery folder:')
        self.l_naip.grid(sticky=E, row=3, column=1, pady=pad)
        self.e_naip = ttk.Entry(root)
        self.e_naip.insert(END, '')
        self.e_naip.grid(row=3, column=2, pady=pad)
        self.b_naip = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_naip, select='folder'))
        self.b_naip.grid(sticky=W, row=3, column=3, pady=pad)

        self.l_ndvi = ttk.Label(root, text='NDVI vegetation threshold:')
        self.l_ndvi.grid(sticky=E, row=4, column=1, pady=pad)
        self.e_ndvi = Entry(root)
        self.e_ndvi.grid(sticky=E, row=4, column=2, pady=pad)
        self.e_ndvi.insert(0, 0.40)
        self.e_ndvi.grid(sticky=E, row=4, column=2, pady=pad)

        self.l_aoi = ttk.Label(root, text='AOI shapefile:')
        self.l_aoi.grid(sticky=E, row=5, column=1, pady=pad)
        self.e_aoi = ttk.Entry(root)
        self.e_aoi.insert(END, '')
        self.e_aoi.grid(row=5, column=2, pady=pad)
        self.b_aoi = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_aoi, select='file',
                                                                            ftypes=[('Shapefile', '*.shp'),
                                                                                    ('All files', '*')]))
        self.b_aoi.grid(sticky=W, row=5, column=3, pady=pad)


        self.spacer1 = ttk.Label(root, text='')
        self.spacer1.grid(sticky=W, row=6, column=1, pady=pad)
        self.prep_run = ttk.Button(root, text='Run',
                                   command=lambda: lidar_prep(self.e_lasbin.get(), self.e_lidardir.get(),
                                                              self.e_in_spatialref.get(), self.e_naip.get(),
                                                              float(self.e_ndvi.get()), self.e_aoi.get()))
        self.prep_run.grid(sticky=E, row=6, column=2)
        root.grid_rowconfigure(17, minsize=80)

        self.instruct = ttk.Label(root, text='     Verify vegetation mask accuracy after running!')
        self.instruct.grid(sticky=EW, row=6, column=3, pady=pad)

        input_ref_shp = self.e_in_spatialref.get()

        # DEM generation
        ######################################################################
        root = self.tabs['DEM generation']  # LiDAR processing to DEM generation widgets

        def dem_generation(lastoolsdir, lidardir, ground_poly, cores, units_code, keep_orig_pts, coarse_step,
                           coarse_bulge, coarse_spike, coarse_down_spike,
                           coarse_offset, fine_step, fine_bulge, fine_spike,
                           fine_down_spike, fine_offset, out_spatialref,
                           dem_resolution, dem_method, in_spatialref=input_ref_shp):
            """This function used LAStools to generate LAS file tiles prepresenitng bare ground, and then converts them
            to a high resolution DEM (1m resolution is default)"""

            # We carry input spatial ref over from the above process, but we should still convert from shp to ref object

            print('DEM generation function running...ADD interpolation optionality')

        self.l_lasbin = ttk.Label(root, text='LAStools /bin/ directory:')
        self.l_lasbin.grid(sticky=E, row=0, column=1)
        self.e_lasbin = ttk.Entry(root)
        self.e_lasbin.insert(END, '')
        self.e_lasbin.grid(row=0, column=2)
        self.b_lasbin = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_lasbin, select='folder'))
        self.b_lasbin.grid(sticky=W, row=0, column=3)

        self.l_lidardir = ttk.Label(root, text='LiDAR data directory:')
        self.l_lidardir.grid(sticky=E, row=1, column=1)
        self.e_lidardir = ttk.Entry(root)
        self.e_lidardir.insert(END, '')
        self.e_lidardir.grid(row=1, column=2)
        self.b_lidardir = ttk.Button(root, text='Browse',
                                     command=lambda: browse(root, self.e_lidardir, select='folder'))
        self.b_lidardir.grid(sticky=W, row=1, column=3)

        self.l_ground_shp = ttk.Label(root, text='Ground area .shp file (optional):')
        self.shp_var = StringVar()
        self.l_ground_shp.grid(sticky=E, row=2, column=1)
        self.e_ground_shp = ttk.Entry(root, textvariable=self.shp_var)
        self.e_ground_shp.grid(row=2, column=2)
        self.b_ground_shp = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_ground_shp,
                                                                                   select='file',
                                                                                   ftypes=[('Shapefile', '*.shp'),
                                                                                           ('All files', '*')]))
        self.b_ground_shp.grid(sticky=W, row=2, column=3)

        self.l_out_spatialref = ttk.Label(root, text='Output spatial reference (.shp)')
        self.l_out_spatialref.grid(sticky=E, row=3, column=1)
        self.e_out_spatialref = ttk.Entry(root)
        self.e_out_spatialref.insert(END, '')
        self.e_out_spatialref.grid(row=3, column=2, pady=pad)
        self.b_out_spatialref = ttk.Button(root, text='Browse',
                                           command=lambda: browse(root, self.e_out_spatialref, select='file',
                                                                  ftypes=[('Shapefile', '*.shp'),
                                                                          ('All files', '*')]))
        self.b_out_spatialref.grid(sticky=W, row=3, column=3, pady=pad)

        # if no ground shapefile is provided, disable the fine setting and just run on "coarse"
        def trace_choice(*args):
            fine_entries = [self.e_f_step, self.e_f_bulge, self.e_f_spike, self.e_f_dspike, self.e_f_offset]
            if self.shp_var.get() == '':
                for widget in fine_entries:
                    widget.config(state=DISABLED)
            else:
                for widget in fine_entries:
                    widget.config(state='normal')

        self.shp_var.trace('w', trace_choice)

        # specify lasground_new parameters

        root.grid_rowconfigure(5, minsize=80)

        self.l_coarse_class = ttk.Label(root, text='standard/coarse classification parameters:')
        self.l_coarse_class.grid(row=5, column=0, columnspan=2)

        self.l_c_step = ttk.Label(root, text='step size:')
        self.l_c_step.grid(sticky=E, row=6)
        self.e_c_step = ttk.Entry(root)
        self.e_c_step.grid(row=6, column=1)

        self.l_c_bulge = ttk.Label(root, text='bulge:')
        self.l_c_bulge.grid(sticky=E, row=7)
        self.e_c_bulge = ttk.Entry(root)
        self.e_c_bulge.grid(row=7, column=1)

        self.l_c_spike = ttk.Label(root, text='spike:')
        self.l_c_spike.grid(sticky=E, row=8)
        self.e_c_spike = ttk.Entry(root)
        self.e_c_spike.grid(row=8, column=1)

        self.l_c_dspike = ttk.Label(root, text='down spike:')
        self.l_c_dspike.grid(sticky=E, row=9)
        self.e_c_dspike = ttk.Entry(root)
        self.e_c_dspike.grid(row=9, column=1)

        self.l_c_offset = ttk.Label(root, text='offset:')
        self.l_c_offset.grid(sticky=E, row=10)
        self.e_c_offset = ttk.Entry(root)
        self.e_c_offset.grid(row=10, column=1)

        self.l_fine_class = ttk.Label(root, text='fine classification parameters (in ground area):')
        self.l_fine_class.grid(row=5, column=2, columnspan=2)

        self.l_f_step = ttk.Label(root, text='step size:')
        self.l_f_step.grid(sticky=E, row=6, column=2)
        self.e_f_step = ttk.Entry(root, state=DISABLED)
        self.e_f_step.grid(row=6, column=3)

        self.l_f_bulge = ttk.Label(root, text='bulge:')
        self.l_f_bulge.grid(sticky=E, row=7, column=2)
        self.e_f_bulge = ttk.Entry(root, state=DISABLED)
        self.e_f_bulge.grid(row=7, column=3)

        self.l_f_spike = ttk.Label(root, text='spike:')
        self.l_f_spike.grid(sticky=E, row=8, column=2)
        self.e_f_spike = ttk.Entry(root, state=DISABLED)
        self.e_f_spike.grid(row=8, column=3)

        self.l_f_dspike = ttk.Label(root, text='down spike:')
        self.l_f_dspike.grid(sticky=E, row=9, column=2)
        self.e_f_dspike = ttk.Entry(root, state=DISABLED)
        self.e_f_dspike.grid(row=9, column=3)

        self.l_f_offset = ttk.Label(root, text='offset:')
        self.l_f_offset.grid(sticky=E, row=10, column=2)
        self.e_f_offset = ttk.Entry(root, state=DISABLED)
        self.e_f_offset.grid(row=10, column=3)

        # specify units
        self.l_lidar_units = ttk.Label(root, text='Units:')
        self.l_lidar_units.grid(sticky=W, row=11, column=2)
        root.grid_rowconfigure(11, minsize=30)
        self.lidar_units = StringVar()
        self.r_lidar_meters = ttk.Radiobutton(root, text='Meters', variable=self.lidar_units, value=' ')
        self.r_lidar_meters.grid(sticky=E, row=12, column=1)
        self.r_lidar_feet = ttk.Radiobutton(root, text='US Feet', variable=self.lidar_units,
                                            value=' -feet -elevation_feet ')
        self.r_lidar_feet.grid(row=12, column=2, pady=pad)
        self.lidar_units.set(' ')

        # specify number of cores
        self.l_lidar_cores = ttk.Label(root, text='Number of cores for processing:')
        self.l_lidar_cores.grid(sticky=E, row=13, column=1, columnspan=2)
        root.grid_rowconfigure(13, minsize=30)
        self.core_num = IntVar()
        self.r1_lidar = ttk.Radiobutton(root, text='1', variable=self.core_num, value=1)
        self.r1_lidar.grid(sticky=E, row=14, column=1)
        self.r2_lidar = ttk.Radiobutton(root, text='2', variable=self.core_num, value=2)
        self.r2_lidar.grid(row=14, column=2)
        self.r4_lidar = ttk.Radiobutton(root, text='4', variable=self.core_num, value=4)
        self.r4_lidar.grid(sticky=W, row=14, column=3)
        self.r8_lidar = ttk.Radiobutton(root, text='8', variable=self.core_num, value=8)
        self.r8_lidar.grid(sticky=E, row=15, column=1)
        self.r16_lidar = ttk.Radiobutton(root, text='16', variable=self.core_num, value=16)
        self.r16_lidar.grid(row=15, column=2)
        self.r32_lidar = ttk.Radiobutton(root, text='32', variable=self.core_num, value=32)
        self.r32_lidar.grid(sticky=W, row=15, column=3, pady=pad)
        self.core_num.set(16)

        self.l_keep_orig_lidar = ttk.Label(root, text='Keep original ground/veg points: ')
        self.l_keep_orig_lidar.grid(sticky=E, row=16, column=1)
        self.keep_orig_lidar = BooleanVar()
        self.c_keep_orig_lidar = ttk.Checkbutton(root, variable=self.keep_orig_lidar)
        self.c_keep_orig_lidar.grid(sticky=W, row=16, column=2)
        self.keep_orig_lidar.set(True)

        units = self.lidar_units.get()[1:-1]

        if units == '':  # Maybe adjust base
            default = 1
        else:
            default = 3.28

        self.l_dem_res = ttk.Label(root, text='DEM resolution (meters):')
        self.l_dem_res.grid(sticky=E, row=17, column=2)
        self.e_dem_res = ttk.Entry(root)
        self.e_dem_res.grid(sticky=W, row=17, column=3, pady=pad)
        self.e_dem_res.insert(0, default)

        self.l_dem_meth = ttk.Label(root, text='Select interpolation method:')
        self.l_dem_meth.grid(sticky=E, row=18, column=2)
        # Linear and natural neighbors refer to TIN based methods, be sure to document
        methods = ['AVERAGE', 'NEAREST', 'SIMPLE', 'Linear', 'Natural neighbors']
        dem_meth = StringVar()
        self.e_dem_meth = ttk.OptionMenu(root, dem_meth, *methods)
        self.e_dem_meth.grid(sticky=W, row=18, column=3)

        # make 'Run' ttk.Button in GUI to call the process_lidar() function
        self.b_lidar_run = ttk.Button(root, text='    Run    ',
                                      command=lambda: dem_generation(lastoolsdir=self.e_lasbin.get(),
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
                                                                     out_spatialref=self.e_out_spatial.get(),
                                                                     dem_resolution=self.e_dem_res.get(),
                                                                     dem_method=dem_meth.get()))

        self.b_lidar_run.grid(sticky=W, row=20, column=2)
        root.grid_rowconfigure(20, minsize=80)

        # Generate thalweg centerline and extract elevation profile
        ######################################################################
        root = self.tabs['Thalweg centerline']

        def detrend_prep(raster_name, flow_polygon, spatial_extent, ft_spatial_ref='', ft_spacing=3, use_filtered_ras=False,
                         centerline_verified=False):
            """Dummy function until we polish up the real one"""
            print('in the function')

        self.remind = ttk.Label(root, text='Create upstream flow polygon in ArcMap')
        self.remind.grid(sticky=E, row=0, column=0)

        self.l_flow_poly = ttk.Label(root, text='Upstream flow shapefile:')
        self.l_flow_poly.grid(sticky=E, row=1, column=0, pady=pad)
        self.e_flow_poly = ttk.Entry(root)
        self.e_flow_poly.grid(sticky=E, row=1, column=1, pady=pad)
        self.e_flow_poly.insert(END, '')
        self.e_flow_poly.grid(row=1, column=1, pady=pad, padx=5)
        self.b_flow_poly = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_flow_poly, select='file',
                                                                                  ftypes=[('Shapefile', '*.shp'), ('All files', '*')]))
        self.b_flow_poly.grid(sticky=W, row=1, column=2, pady=pad)

        self.l_extent = ttk.Label(root, text='Extent shapefile:')
        self.l_extent.grid(sticky=E, row=2, column=0, pady=pad)
        self.e_extent = ttk.Entry(root)
        self.e_extent.grid(row=2, column=1, pady=pad)
        self.e_extent.insert(END, '')
        self.e_extent.grid(row=2, column=1, pady=pad, padx=5)
        self.b_extent = ttk.Button(root, text='Browse', command=lambda: browse(root, self.e_extent, select='file',
                                                                            ftypes=[('Shapefile', '*.shp'),
                                                                                    ('All files', '*')]))
        self.b_extent.grid(sticky=W, row=2, column=2, pady=pad)

        self.l_filt = ttk.Label(root, text='Filter passes (15x default):')
        self.l_filt.grid(sticky=E, row=3, column=0, pady=pad)
        self.e_filt = Entry(root)
        self.e_filt.grid(sticky=E, row=3, column=1, pady=pad)
        self.e_filt.insert(END, 15)
        self.e_filt.grid(row=3, column=1, pady=pad, padx=5)

        self.l_smooth = ttk.Label(root, text='Smoothing distance (DEM units):')
        self.l_smooth.grid(sticky=E, row=4, column=0, pady=pad)
        self.e_smooth = Entry(root)
        self.e_smooth.grid(sticky=E, row=4, column=1, pady=pad)
        self.e_smooth.insert(END, 20)
        self.e_smooth.grid(row=4, column=1, pady=pad, padx=5)

        # DUMMY VARIABLES IN FUNCTION JUST WORKING ON LAYOUT RN
        self.b_detrend_prep1 = ttk.Button(root, text='    Run    ',
                                      command=lambda: detrend_prep(raster_name, flow_polygon=self.e_flow_poly.get(),
                                                                   spatial_extent=self.e_extent.get(), ft_spatial_ref='', ft_spacing=3,
                                                                   use_filtered_ras=False, centerline_verified=False))



if __name__ == '__main__':
    gcs_gui().mainloop()
