import arcpy
import os
import pandas
import pandas as pd
from typing import List, Tuple, Union
from matplotlib import pyplot as plt
import numpy as np
from file_functions import err_info, spatial_license

# Define detrending functions
######################################################################


def prep_xl_file(
    xyz_csv: str,
    in_columns: List[str] = ['LOCATION', 'POINT_X', 'POINT_Y', 'Value'],
) -> Tuple[np.array, np.array, str]:
    """"This function takes the .csv file exported during detrending prep and returns arrays representing the longitudinal
    thalweg elevation profile.
    Returns: (List) [array of distance downstream, array of thalweg z values, the original xyz csv file]"""
    list_of_lists = [None, None, None, None]

    elevation_df = pd.read_csv(xyz_csv)
    for j, header in enumerate(in_columns):
        list_of_lists[j] = elevation_df.loc[:, [header]].squeeze().to_numpy()

    location = np.int_(list_of_lists[0])
    z = np.around(list_of_lists[-1], 9)

    return (location, z, xyz_csv)


def linear_fit(
    location_np: np.array,
    z_np: np.array,
    xyz_table_loc: str,
    bp_list: List[Union[int, float, None]] = [],
) -> Tuple[List[List[float]], np.array, np.array, float]:
    """Applies a linear fit to piecewise sections of the longitudinal profile, each piece is stored in split_list"""

    print("Applying linear fit...")

    # Initiate lists for breakpoint based linear fit
    split_locs_list = []
    split_z_list = []
    fit_params = []
    z_fit_list = []

    if len(bp_list) != 0:
        bp_list.insert(0, 0)
        bp_list.append(int(location_np[-1]))
        print("Breakpoints imported...")
    else:
        print("No breakpoint imported...")

    # Set up arrays and calculate point spacing
    location_np = np.int_(location_np)
    point_spacing = int(location_np[1]) - int(location_np[0])
    z_np = np.float_(z_np)
    z_np = np.around(z_np, 9)  # Round z to 9 decimal points

    if len(bp_list) > 0:
        bp_indices = [int(dist / point_spacing) for dist in bp_list]

        for i in bp_indices[1:]:
            index = bp_indices.index(i)

            if i == 1:
                split_locs_list.append(location_np[:i])
                split_z_list.append(z_np[:i])

            elif i != bp_indices[-1]:
                split_locs_list.append(location_np[bp_indices[index - 1]:i])
                split_z_list.append(z_np[bp_indices[index - 1]:i])

            elif i == bp_indices[-1]:
                split_locs_list.append(location_np[bp_indices[index - 1]:])
                split_z_list.append(z_np[bp_indices[index - 1]:])
        print("Breakpoints added")

        # Get fit parameters for each section of the data
        if len(split_z_list) == len(split_locs_list):
            for i in range(len(split_locs_list)):
                temp_locs = np.int_(np.array(split_locs_list[i]))
                temp_zs = np.float_(np.array(split_z_list[i]))

                m, b = np.polyfit(temp_locs, temp_zs, 1)
                fit_params.append([m, b])

        else:
            print("Something went wrong, list lengths do not match...")

        # Make list storing the lengths of each breakpoint split location/z array segment
        lengths = []
        for i in range(len(split_z_list)):
            lengths.append(len(split_z_list[i]))

        # For each segment, generate an array storing the fitted z values
        i = 0
        while i < len(lengths):
            for j in range(lengths[i]):
                z_fit_list.append(
                    split_locs_list[i][j] * fit_params[i][0] + fit_params[i][1])
            i += 1

    else:
        m, b = np.polyfit(location_np, z_np, 1)
        fit_params = [[m, b]]
        print("Fit params [m, b]: " + str(fit_params))

        location_list = location_np.tolist()
        for j in location_list:
            z_fit_list.append(j*fit_params[0][0]+fit_params[0][1])

    # Save fitted values to the output .csv table
    elevation_df = pd.read_csv(xyz_table_loc)
    elevation_df['z_fit'] = np.array(z_fit_list)
    elevation_df.to_csv(xyz_table_loc)

    # Calculate residual
    residual = []
    for i in range(len(z_fit_list)):
        residual.append(z_np[i] - z_fit_list[i])
    mean_z = (sum(z_np) / len(z_np))
    squared_real = []
    squared_res = []

    # Calculate the R^2 value
    for i in range(len(residual)):
        squared_real.append((z_np[i] - mean_z) ** 2)
        squared_res.append(residual[i] ** 2)

    r_squared = 1 - (sum(squared_res) / sum(squared_real))

    # Convert residual and z fit values to an array
    z_fit = np.array(z_fit_list)
    residual = np.array(residual)
    print('Done')

    return (fit_params, z_fit, residual, r_squared)

@err_info
@spatial_license
def detrend_that_raster(
    xyz_csv: str,
    in_dem: str,
    aoi_shp: str = '',
) -> str:
    """Generates a detrended DEM from a the fitted xyz .csv file and an input .tif dem"""
    # Set up directory structure and environment
    out_dir = os.path.dirname(xyz_csv)
    temp_files = out_dir + '\\temp_files'
    if not os.path.exists(temp_files):
        os.makedirs(temp_files)

    arcpy.env.workspace = temp_files
    arcpy.overwriteoutput = True

    out_dem = out_dir + '\\ras_detren.tif'
    spatial_ref = arcpy.Describe(in_dem).spatialReference
    arcpy.env.extent = arcpy.Describe(in_dem).extent

    # Create dataframe storing the fitted xyz .csv values
    fit_col = 'z_fit'
    cols = ['FID', 'Shape', 'POINT_X', 'POINT_Y', fit_col]
    xyz_df = pandas.read_csv(xyz_csv, usecols=cols)

    # Generate station points with fitted z values
    points = arcpy.MakeXYEventLayer_management(
        xyz_csv,
        "POINT_X",
        "POINT_Y",
        out_layer='fit_station_points',
        spatial_reference=spatial_ref,
        in_z_field=fit_col,
    )
    points = arcpy.SaveToLayerFile_management(points, 'fit_station_points.lyr')
    points = arcpy.CopyFeatures_management(points)

    # Delete non-relevant csv columns
    fields = [f.name for f in arcpy.ListFields(points)]
    fields2delete = list(set(fields) - set(cols))
    points = arcpy.DeleteField_management(points, fields2delete)

    # Calculate dem cell size and generate thiessen raster from fitted station points
    print("Creating Thiessen polygons...")
    cell_size1 = arcpy.GetRasterProperties_management(in_dem, "CELLSIZEX")
    cell_size = float(cell_size1.getOutput(0))

    thiessen = arcpy.CreateThiessenPolygons_analysis(
        points,
        "thiespoly.shp",
        fields_to_copy='ALL',
    )

    z_fit_ras = arcpy.PolygonToRaster_conversion(
        thiessen,
        fit_col,
        'theis_ras.tif',
        cell_assignment="MAXIMUM_AREA",
        cellsize=cell_size,
    )

    # Detrend in_dem by subtracting thiessen raster values from it
    detrended_dem = arcpy.Raster(in_dem) - arcpy.Raster(z_fit_ras)

    if aoi_shp == '':
        detrended_dem.save(out_dem)
    else:
        no_clip = temp_files + '\\ras_dt_nc.tif'
        detrended_dem.save(no_clip)
        arcpy.Clip_management(
            no_clip,
            out_raster=out_dem,
            in_template_dataset=aoi_shp,
            clipping_geometry='ClippingGeometry',
        )

    return out_dem


# Define plotting functions
######################################################################

def diagnostic_quick_plot(
    location_np: np.array,
    z_np: np.array,
    out_dir: str,
) -> str:
    """Generates a basic plot showing the thalweg elevation profile.
    A folder where plots can be saved (sub-folder generated)."""
    x_plot = location_np
    y_plot = z_np
    plt.plot(x_plot, y_plot, 'r', label='Thalweg elevation profile')

    # Define plotting extent
    plt.xlim(min(location_np), max(location_np))
    plt.ylim(min(z_np), max(z_np))

    # Set up plot labels
    plt.xlabel('Thalweg distance downstream', fontsize='small')
    plt.ylabel('Thalweg elevation', fontsize='small')

    # Format plot
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.minorticks_on()
    plt.xticks(fontsize='x-small')
    plt.yticks(fontsize='x-small')

    # Save plot, return address
    fig = plt.gcf()
    fig.set_size_inches(6, 3)

    out_png = out_dir + '\\thalweg_z_plot.png'
    plt.savefig(
        out_png,
        dpi=300,
        bbox_inches='tight',
    )
    plt.cla()

    return out_png


def linear_fit_plot(
    location_np: np.array,
    z_np: np.array,
    fit_params: List[List[float]],
    fit_np: np.array,
    out_dir: str,
) -> str:
    """Generates a plot showing linear fit models between breakpoints
    Inputs: Numpy array of distance downstream, thalweg z values. List of fit parameters output from detrending function.
    A folder where plots can be saved (sub-folder generated)."""

    # Prep input arrays
    y1_plot = z_np
    y2_plots = []

    for sub in fit_params:
        y2_plots.append(sub[0]*location_np + sub[1])

    # Initiate plot
    plt.plot(location_np, y1_plot, 'r', label="Thalweg elevation profile")
    plt.plot(location_np, fit_np, 'b',
             label='Piecewise linear fit', linewidth=0.75)

    # Define plotting extent
    plt.xlim(min(location_np), max(location_np))
    plt.ylim(min(z_np), max(z_np))

    # Set up plot labels
    plt.xlabel('Thalweg distance downstream', fontsize='small')
    plt.ylabel('Thalweg elevation', fontsize='small')
    plt.title('Linear piecewise fit', fontsize='small')

    # Format plot
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.minorticks_on()
    plt.xticks(fontsize='x-small')
    plt.yticks(fontsize='x-small')
    plt.legend(loc=1, fontsize='x-small')

    # Save plot, return address
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    out_png = out_dir + '\\fit_plot.png'
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.cla()

    return out_png


def make_residual_plot(
    location_np: np.array,
    residual_np: np.array,
    r2: float,
    out_dir: str,
) -> str:
    """Plots residuals across the longitudinal profile, shows the R^2 value. Outlier values removed for view-ability.
    Inputs: Numpy arrays of distance downstream and fit residuals. A R-squared value (float).
    A folder where plots can be saved (sub-folder generated).
    Returns: A figure showing linear fit residuals."""

    # Prep input arrays and initiate plot
    y_zero = 0*location_np
    plt.scatter(location_np, residual_np, s=1, c='r')
    plt.plot(location_np, y_zero, c='b')

    # Define plotting extent using percentile (removes outliers for improved plot  output)
    bottom = np.percentile(residual_np, 1)
    top = np.percentile(residual_np, 99)
    plt.xlim(0, max(location_np))
    plt.ylim(bottom, top)

    # Set up plot labels
    plt.xlabel('Thalweg distance downstream', fontsize='small')
    plt.ylabel('Fit residual', fontsize='small')
    plt.title('Residuals: R^2 = %.4f' % r2, fontsize='small')

    # Format plot
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.minorticks_on()
    plt.xticks(fontsize='x-small')
    plt.yticks(fontsize='x-small')

    # Save plot, return address
    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    out_png = out_dir + '\\residual_plot.png'
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.cla()

    return out_png


def fit_params_txt(
    fit_params: List[List[float]],
    bp_list: List[Union[int, float]],
    out_dir: str,
) -> str:
    """Generates a text file in the same folder as the detrending plots that lists applied linear fit equations"""

    # Create .txt file and copy breakpoint list
    text_dir = out_dir + '\\detrending_fit_eqs.txt'
    text_file = open(text_dir, 'w+')
    bps_form = [i for i in bp_list]

    # Write to and save .txt file
    for count, params in enumerate(fit_params):
        if len(bp_list) != 0:
            text_file.write('From %s to %s: %.4f * dist_downstream + %.4f\n' %
                            (bps_form[count], bps_form[count+1], params[0], params[1]))
        else:
            text_file.write(
                'For full reach: %.4f * dist_downstream + %.4f\n' % (params[0], params[1]))
    text_file.close()

    return text_dir
