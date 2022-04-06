import csv
import os
import math
import numpy as np
import pandas as pd
import scipy
from scipy.stats import variation
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import matplotlib.colors as colors_module
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.cm
import plotly
import plotly.graph_objects as go
import plotly.express as pex
from itertools import combinations
import seaborn as sns

import file_functions
import openpyxl as xl


def gcs_plotter(detrended_dem, zs, fields=['Ws', 'Zs', 'Ws_Zs'], together=False):
    """This function makes longitudinal profile plots for given fields across each key z saving them to a folder.
     If aligned_table is defined as the aligned csv, plots showing each key z profile as sub-plots for a given field are saved as well."""

    if detrended_dem == '':
        print('Error: Must input detrended DEM parameter in the GUI to set up output folder location')
        return

    if type(zs) == str:
        zs = file_functions.string_to_list(zs, format='float')
    elif type(zs) != list:
        print(
            'Error: Key flow stage parameter input incorrectly. Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
    if len(zs) == 0:
        print('Error: Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
        return
    zs.sort()

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    gcs_dir = dem_dir + '\\gcs_tables'

    # get units for labeling
    u = file_functions.get_label_units(detrended_dem)[0]

    # define output directory
    if together:
        out_dir = dem_dir + '\\nesting_analysis'
    else:
        out_dir = dem_dir + '\\stage_analysis'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    colors = ['black', 'blue', 'grey', 'orange', 'red']
    landforms = ['Oversized', 'Const. Pool', 'Normal', 'Wide Bar', 'Nozzle']

    # create subplots for each flow stage for a given gcs series (i.e. Ws, Zs, Ws*Zs)
    if together:
        for field in fields:
            xs = []
            ys = []  # previosuly was [] for i in zs
            full_ys = []
            for count, z in enumerate(zs):
                ys.append([])
                label = file_functions.float_keyz_format(z) + u
                table_loc = gcs_dir + '\\%s_gcs_table.csv' % label
                table_df = pd.read_csv(table_loc)

                table_df.sort_values('dist_down', inplace=True)
                xs.append(table_df.loc[:, 'dist_down'].to_numpy())
                codes = table_df.loc[:, 'code'].to_numpy()

                for num in range(-2, 3):  # Make arrays representing each landform type
                    y_temp = table_df.loc[:, field].to_numpy()
                    y_temp[codes != num] = np.nan
                    ys[count].append(y_temp)
                    table_df = pd.read_csv(table_loc)
                    table_df.sort_values('dist_down', inplace=True)

                full_ys.append(table_df.loc[:, field].to_numpy())

            fig, ax = plt.subplots(len(zs), sharey=True)
            fig.subplots_adjust(hspace=0.4)
            fig_name = out_dir + '\\%s_nesting_gcs_plots.png' % field
            ax[0].set_title('%s series' % field)

            for count, z in enumerate(zs):
                x = xs[count]
                ymax = 0
                ax[count].plot(x, full_ys[count], color=colors[2])
                for i, y in enumerate(ys[count]):
                    ax[count].plot(x, y, color=colors[i], label=landforms[i])
                    temp_max = np.amax(np.array([np.abs(np.nanmin(y)), np.abs(np.nanmax(y))]))
                    if temp_max >= ymax and ymax <= 5:
                        ymax = math.ceil(temp_max)
                    elif ymax > 5:
                        ymax = 5
                ax[count].set_ylim(-1 * ymax, ymax)
                ax[count].set_ylabel(field)
                ax[count].set_yticks(np.arange(-1 * ymax, ymax, 1), minor=False)
                ax[count].grid(True, which='both', color='gainsboro', linestyle='--')
                ax[count].set_xlim(0.0, np.max(x))
                ax[count].set_xticks(np.arange(0, np.max(x), 250))

            ax[count].set_xlabel('Thalweg distance downstream (ft)')
            ax[count].legend(loc='lower center', ncol=len(landforms), fontsize=8)  # Adds legend to the bottom plot
            fig.set_size_inches(12, 6)
            plt.savefig(fig_name, dpi=300, bbox_inches='tight')
            plt.cla()
        plt.close('all')

    # create subplots for each gcs series (i.e. Ws, Zs, Ws*Zs) for a given flow stage
    elif not together:
        for z in zs:
            # get data for the flow stage
            label = file_functions.float_keyz_format(z) + u
            table_loc = gcs_dir + '\\%s_gcs_table.csv' % label
            table_df = pd.read_csv(table_loc)

            # extract distance downstream for the x axis of each subplot, as well as landform codes
            table_df.sort_values('dist_down', inplace=True)
            x = table_df.loc[:, 'dist_down'].to_numpy()
            codes = table_df.loc[:, 'code'].to_numpy()

            ys = [[] for i in fields]  # previosuly was [] for i in zs
            full_ys = []
            for count, field in enumerate(fields):
                for num in range(-2, 3):  # Make arrays representing each landform type
                    y_temp = table_df.loc[:, field].to_numpy()
                    y_temp[codes != num] = np.nan
                    ys[count].append(y_temp)
                    table_df = pd.read_csv(table_loc)
                    table_df.sort_values('dist_down', inplace=True)

                full_ys.append(table_df.loc[:, field].to_numpy())

            fig, ax = plt.subplots(len(fields), sharey=True)
            fig.subplots_adjust(hspace=0.4)
            fig_name = out_dir + '\\%s_gcs_plots.png' % label
            ax[0].set_title('%s stage' % label)

            for count, field in enumerate(fields):
                ymax = 0
                ax[count].plot(x, full_ys[count], color=colors[2])
                for i, y in enumerate(ys[count]):
                    ax[count].plot(x, y, color=colors[i], label=landforms[i])
                    temp_max = np.amax(np.array([np.abs(np.nanmin(y)), np.abs(np.nanmax(y))]))
                    if temp_max >= ymax and ymax <= 5:
                        ymax = math.ceil(temp_max)
                    elif ymax > 5:
                        ymax = 5

                ax[count].set_ylim(-1 * ymax, ymax)
                ax[count].set_ylabel(field)
                ax[count].set_yticks(np.arange(-1 * ymax, ymax, 1), minor=False)
                ax[count].grid(True, which='both', color='gainsboro', linestyle='--')
                ax[count].set_xlim(0.0, np.max(x))
                ax[count].set_xticks(np.arange(0, np.max(x), 250))

            ax[count].set_xlabel('Thalweg distance downstream (ft)')
            ax[count].legend(loc='lower center', ncol=len(landforms), fontsize=8)  # Adds legend to the bottom plot
            fig.set_size_inches(12, 6)
            plt.savefig(fig_name, dpi=300, bbox_inches='tight')
            plt.cla()
        plt.close('all')
    print('GCS plots saved @ %s' % out_dir)


def heat_plotter(detrended_dem, zs, together=False):
    """IN: Detrended DEM path (str), a string or list with flow stage height floats, whether to plot stages together or separate (boolean).
    Returns: png heat-plot(s) for each input flow stage either all on one plot (together=True), or separately (together=False)"""

    if detrended_dem == '':
        print('Error: Must input detrended DEM parameter in the GUI to set up output folder location')
        return

    if type(zs) == str:
        zs = file_functions.string_to_list(zs, format='float')
    elif type(zs) != list:
        print(
            'Error: Key flow stage parameter input incorrectly. Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
    if len(zs) == 0:
        print('Error: Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
        return

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    gcs_dir = dem_dir + '\\gcs_tables'

    # get units for labeling
    u = file_functions.get_label_units(detrended_dem)[0]

    # use [[zs]] or [[z1], [z2]] structure to control plotting
    if together:
        top_zs = [zs]
        out_dir = dem_dir + '\\nesting_analysis'
        title = out_dir + '\\stages_heatplots.png'
    else:
        top_zs = []
        for i in zs:
            top_zs.append([i])
        out_dir = dem_dir + '\\stage_analysis'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # set up subplot arrangement
    for sub_zs in top_zs:
        if len(sub_zs) > 6:
            print('Warning: > 6 flow stages input for nesting analysis, may seriously impair plot/analysis quality!')
            fig, axs = plt.subplots(ncols=int(len(sub_zs)), figsize=(10, 3))
        elif len(sub_zs) == 4 or len(sub_zs) == 6:
            fig, axs = plt.subplots(nrows=2, ncols=int(len(sub_zs) / 2), figsize=(10, 3))
        else:
            fig, axs = plt.subplots(ncols=int(len(sub_zs)), figsize=(10, 3))

        fig.subplots_adjust(hspace=0.5, wspace=0.3, left=0.07, right=0.93)

        for count, ax in enumerate(axs):
            z = sub_zs[count]
            label = file_functions.float_keyz_format(z) + u

            if not together:
                title = out_dir + '\\%s_heatplot.png' % label

            # create heat-plots
            data = pd.read_csv(gcs_dir + '\\%s_gcs_table.csv' % label)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, ['Ws']].to_numpy()
            y = data.loc[:, ['Zs']].to_numpy()

            ax.set_aspect('equal', adjustable='box')
            ax.hexbin(x, y, gridsize=30, cmap='YlOrRd', extent=(-3, 3, -3, 3))
            ax.set(xlim=(-3, 3), ylim=(-3, 3))
            ax.axhline(y=0.5, xmin=0, xmax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axhline(y=0.5, xmin=0.583, xmax=1, color='#9e9e9e', linestyle='--')
            ax.axhline(y=-0.5, xmin=0, xmax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axhline(y=-0.5, xmin=0.583, xmax=1, color='#9e9e9e', linestyle='--')

            ax.axvline(x=-0.5, ymin=0, ymax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axvline(x=-0.5, ymin=0.583, ymax=1, color='#9e9e9e', linestyle='--')
            ax.axvline(x=0.5, ymin=0, ymax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axvline(x=0.5, ymin=0.583, ymax=1, color='#9e9e9e', linestyle='--')

            ax.text(0.20, 0.05, 'Const. Pool', horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes)
            ax.text(0.18, 0.95, 'Nozzle', horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes)
            ax.text(0.82, 0.95, 'Wide Bar', horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes)
            ax.text(0.82, 0.05, 'Oversized', horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes)

            ax_title = '%s Ws, Zs heat-plot' % label
            ax.set_title(ax_title)
            ax.set_xlabel('Standardized width (Ws)')
            ax.set_ylabel('Standardized detrended elevation (Zs)')

        # save .png
        plt.savefig(title, dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        print('A plot with heat-plots for all stage heights %s is @ %s' % (zs, title))

    return out_dir


def landform_pie_charts(detrended_dem, zs, together=False):
    """"""
    labels = ['Oversized', 'Const.Pool', 'Normal', 'Wide Bar', 'Nozzle']
    colors = ['black', 'blue', 'grey', 'orange', 'red']

    if detrended_dem == '':
        print('Error: Must input detrended DEM parameter in the GUI to set up output folder location')
        return

    if type(zs) == str:
        zs = file_functions.string_to_list(zs, format='float')
    elif type(zs) != list:
        print(
            'Error: Key flow stage parameter input incorrectly. Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
    if len(zs) == 0:
        print('Error: Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')
        return

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    gcs_dir = dem_dir + '\\gcs_tables'

    # get units for labeling
    u = file_functions.get_label_units(detrended_dem)[0]

    # use [[z1,z2,z3]] or [[z1], [z2]] structure to control plotting
    if together:
        top_zs = [zs]
        out_dir = dem_dir + '\\nesting_analysis'
        title = out_dir + '\\stages_landform_pies.png'
    else:
        top_zs = []
        for i in zs:
            top_zs.append([i])
        out_dir = dem_dir + '\\stage_analysis'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # set up subplot arrangement
    for sub_zs in top_zs:
        if len(sub_zs) > 6:
            print('Warning: > 6 flow stages input for nesting analysis, may seriously impair plot/analysis quality!')
            fig, axs = plt.subplots(ncols=int(len(sub_zs)), figsize=(10, 3))
        elif len(sub_zs) == 4 or len(sub_zs) == 6:
            fig, axs = plt.subplots(nrows=2, ncols=int(len(sub_zs) / 2), figsize=(10, 3))
        else:
            fig, axs = plt.subplots(ncols=int(len(sub_zs)), figsize=(10, 3))

        fig.subplots_adjust(hspace=0.5, wspace=0.3, left=0.07, right=0.93)

        percents = []
        middle_index = math.trunc(len(sub_zs)/2)

        for count, ax in enumerate(axs):
            z = sub_zs[count]
            label = file_functions.float_keyz_format(z) + u

            # record total occurrences for each land form code [-2, -1, 0, 1, 2] and calculate percents
            z_df = pd.read_csv(gcs_dir + '\\%s_gcs_table.csv' % label)
            codes = z_df.loc[:, 'code'].to_numpy()
            total = len(codes)
            counts = []
            temp_percents = []
            for num in range(-2, 3):
                count = (codes == num).sum()
                counts.append(count)
                temp_percents.append((count / total) * 100)
                percents.append(np.array(temp_percents))

            ax.pie(percents[count], labels=labels, labeldistance=None, autopct='%1.1f%%', textprops={'color': "w"}, colors=colors)
            ax.set_title(label)
            ax.title.set_position([0.5, 0.92])

        # Adds legend to the bottom middle of the plot
        axs[middle_index].legend(bbox_to_anchor=(0.32, 0.07), bbox_transform=ax.transAxes, ncol=len(labels), fontsize=8)

        if not together:
            title = out_dir + '\\%s_landform_pies.png' % label

        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        plt.savefig(title, dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        print('Pie plots for landform abundances @ %s' % out_dir)

    return out_dir


def nested_landform_sankey(detrended_dem, zs=[], ignore_normal=False):
    """Creates Sankey diagrams showing nested landform relationships. Can be done across a class, transition occurences are normalized as a % for each reach."""

    if detrended_dem == '':
        print('Error: Must input detrended DEM parameter in the GUI to set up output folder location')
        return

    if type(zs) == str:
        zs = file_functions.string_to_list(zs, format='float')
    elif type(zs) != list:
        print(
            'Error: Key flow stage parameter input incorrectly. Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)')

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    gcs_dir = dem_dir + '\\gcs_tables'
    out_dir = dem_dir + '\\nesting_analysis'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # get units for labeling
    u = file_functions.get_label_units(detrended_dem)[0]

    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}  # code number and corresponding MU
    print('Sankey landform diagram plotting comencing...')

    source = []
    target = []
    value = []

    zs.sort()
    aligned_df = pd.read_csv(gcs_dir + '\\aligned_gcs_table.csv')
    data = aligned_df.dropna()

    # create a list that stores the landform code values for each aligned flow stage series
    code_df_list = []
    for z in zs:
        label = file_functions.float_keyz_format(z) + u

        code_df_temp = data.loc[:, [('code_%s' % label)]].squeeze()
        code_df_list.append(code_df_temp.values.tolist())

    # make lists of tuples storing each side of each cross-sections step wise transition
    transitions = []
    for i, t in enumerate(code_df_list[:-1]):
        transitions.append(list(zip(t, code_df_list[i + 1])))

    # initialize list of lists to count abundance
    unique_nests = [list(set(i)) for i in transitions]
    unique_nest_counts = [list(np.zeros(len(i), dtype=int)) for i in unique_nests]

    # Calculates totals of occurences for each incrementing flow stage transition, ex: 0.2->0.7, 0.7->2.6, 2.6->5.2
    for j, zipped in enumerate(transitions):
        for pair in zipped:
            i = unique_nests[j].index(pair)
            unique_nest_counts[j][i] += 1

    nest_abundances = [list(zip(unique_nests[i], unique_nest_counts[i])) for i in range(len(unique_nests))]

    # set up lists to store Sankey diagram node information
    label_list = []
    x_list = []
    y_list = []
    colors_list = []
    x_num = 0.1

    for z in zs:
        if not ignore_normal:
            label_list.extend(['Oversized', 'Const.Pool', 'Normal', 'Wide Bar', 'Nozzle'])
            x_list.extend([x_num for j in range(5)])
            x_num += 0.4
            y_list.extend([0.1, 0.3, 0.5, 0.7, 0.9])
            colors_list.extend(['black', 'blue', 'grey', 'orange', 'red'])

        else:
            label_list.extend(['Oversized', 'Const.Pool', 'Wide Bar', 'Nozzle'])
            x_list.extend([x_num for j in range(4)])
            x_num += 0.4
            y_list.extend([0.2, 0.4, 0.6, 0.8])
            colors_list.extend(['black', 'blue', 'orange', 'red'])

    nodes = {
        "label": label_list,
        "x": x_list,
        "y": y_list,
        "color": colors_list,
        'pad': 15}

    for j, nests in enumerate(nest_abundances):
        for i in nests:
            if ignore_normal == False:
                index_adjust = 2 + j * 5
                source.append(int(i[0][0] + index_adjust))
                target.append(int(i[0][1] + index_adjust + 5))
                value.append(float(i[1] + index_adjust))
            elif 0 not in i[0]:
                index_adjust = 2 + j * 4
                if int(i[0][0]) < 0:
                    source.append(int(i[0][0] + index_adjust))
                else:
                    source.append(int(i[0][0] + index_adjust - 1))
                if int(i[0][1]) < 0:
                    target.append(int(i[0][1] + index_adjust + 4))
                else:
                    target.append(int(i[0][1] + index_adjust + 3))
                value.append(float(i[1] + index_adjust))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=nodes,  # 10 Pixels
        link={
            "source": source,
            "target": target,
            "value": value}))

    fig.write_image(out_dir + '\\landform_transitions.png', scale=5)
    print('Sankey landform transition plots saved @ %s' % out_dir)



