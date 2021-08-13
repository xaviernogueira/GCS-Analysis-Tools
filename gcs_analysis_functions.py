import arcpy
from arcpy import da
import file_functions
from file_functions import *
import pandas as pd
import os

arcpy.env.overwriteOutput = True


@err_info
# Define functions that create standardized Ws, Zs, and C(Ws, Zs) values along with landform classification codes
##################################################################################################################




def clean_in_table(table, w_field='W', z_field='Z', dist_field='dist_down'):
    """Renames columns corresponding to W, Z, and dist_down, if they are not already columns"""
    check_use(table)
    df = pd.read_csv(table)

    for old_name, replacement_name in [(w_field, 'W'), (z_field, 'Z'), (dist_field, 'dist_down')]:
        if replacement_name not in df.columns.tolist():
            if old_name in df.columns.tolist():
                df.rename(columns={old_name: replacement_name})
                logging.info('Renamed %s to %s in %s' % (old_name, replacement_name, table))
            else:
                logging.exception('Cannot find column named %s or %s in %s' % (old_name, replacement_name, table))

    df.to_csv(table, index=False)

    return df


def standardize(table, fields):
    """Makes standardized version of field in csv table by subtracting each value by mean and dividing by standard deviation.
    Creates a new column Ws_Zs storing C(Ws,Zs) values"""
    check_use(table)
    df = pd.read_csv(table)

    if type(fields) == list:
        for f in fields:
            new_field = fields + 's'
            df[new_field] = (df[f] - np.mean(df[f])) * 1.0 / np.std(df[f])

    df['%s_%s' % (fields[0], fields[1])] = df[fields[0]] * df[fields[1]]
    df.to_csv(table, index=False)

    return df


def landforms(table, zs_field='Zs', ws_field='Ws', na=-9999):
    '''Classifies each row by corresponding landform type:
        oversized, nozzle, constricted pool, wide bar, normal channel
        Adds columns to input table'''

    check_use(table)
    df = pd.read_csv(table)

    df['normal'] = [zs * ws if abs(zs) <= 0.5 or abs(ws) <= 0.5 else na for zs, ws in zip(df[zs_field], df[ws_field])]
    df['wide_bar'] = [zs * ws if (zs > 0.5 and ws > 0.5) else na for zs, ws in
                      zip(df[zs_field], df[ws_field])]
    df['const_pool'] = [zs * ws if (zs < -0.5 and ws < -0.5) else na for zs, ws in
                        zip(df[zs_field], df[ws_field])]
    df['nozzle'] = [zs * ws if (zs > 0.5 and ws < -0.5) else na for zs, ws in
                    zip(df[zs_field], df[ws_field])]
    df['oversized'] = [zs * ws if (zs < 0.5 and ws > 0.5) else na for zs, ws in
                       zip(df[zs_field], df[ws_field])]

    df['code'] = [-2 if df['oversized'][i] != na
                  else -1 if df['const_pool'][i] != na
    else 0 if df['normal'][i] != na
    else 1 if df['wide_bar'][i] != na
    else 2 if df['nozzle'][i] != na
    else 0
                  # Was na, but since for whatever reason normal channel is not mutually exclusive, we are going to hard code this as 0
                  for i in range(len(df))
                  ]

    df["code"].fillna(0, inplace=True)
    df.to_csv(table, index=False)
    return df


def main_classify_landforms(tables, w_field, z_field, dist_field, out_folder):
    """Classifies river segments as normal, wide bar, constricted pool, oversized, or nozzle

    Args:
        tables: a list of attribute table filenames for each set of wetted polygon rectangular XS's
        w_field: name of the attribute table field corresponding to width
        z_field: name of the attribute table field corresponding to detrended bed elevation
        dist_field: name of the attribute table field corresponding to distance downstream

    Returns:
        For each input table:
            a .csv containing dist_down, W, Z, W_s, Z_s, W_s*Z_s, and landform classification/code fields
            adds these computed values to attribute tables of corresponding wetted polygon rectangular XS's
    """
    logging.info('Classifying landforms...')
    out_polys = []
    fields = [w_field, z_field]
    std_fields = [field + '_s' for field in fields]
    std_pairs = list(itertools.combinations(std_fields, 2))

    width_dir = out_folder + "\\analysis_shapefiles"
    list_of_files_width_folder = [f for f in os.listdir(width_dir) if
                                  os.path.isfile(os.path.join(width_dir, f))]  # add width rectangles to a list
    width_files = []
    for file in list_of_files_width_folder:
        if file[:16] == "width_rectangles" and file[-4:] == ".shp" and file[-8:-6] != '_0':
            width_files.append(file)

    for i in range(len(tables)):
        table = tables[i]
        clean_in_table(table, w_field=w_field, z_field=z_field, dist_field=dist_field)
        standardize(table, fields=fields)
        df = landforms(table)
        # df = landform_polygons(table, width_dir + "\\" + width_files[i])
        # table.replace('.csv', '.shp')
        out_polys.append(width_dir + "\\" + width_files[i])

    logging.info('Finished.')

    return out_polys


# Main function that conducts GCS analysis
############################################


def key_zs_gcs(detrended_dem, zs, xs_lengths, spacing, clip_poly=''):
    '''This function does a full GCS analysis using three specific key Zs that can include any float. Results saved
    to the gcs_ready_tables, as well as plotted. Results are aligned to the existing csv to facilitate landform analysis
    detrend
    wetted_folder is the folder containing small increment wetted polygons
    If key_zs parameter is an empty list, a range from 0 to max_stage (deafult is 20) makes gcs csvs at 1ft increments.
    wetted_folder parameter (optional) allows for a specific folder containing wetted polygons to be used instead of the assumed file structure.'''

    # Convert stage heights to list format
    if type(zs) == str:
        zs = string_to_list(zs, format='float')
    if type(xs_lengths) == str:
        xs_lengths = string_to_list(xs_lengths, format='float')
    if type(spacing) == str:
        try:
            spacing = int(spacing)
        except TypeError:
            print('Error: Cross-section spacing parameter must be a integer!')

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)

    if len(dem_dir) == 0:
        print('Error: Please select valid detrended DEM file')
        return

    del_files = []  # stored deletable files

    lines_dir = dem_dir + '\\centerlines'
    wetted_dir = dem_dir + '\\wetted_polygons'
    temp_files = dem_dir + '\\temp_files'
    out_dir = dem_dir + '\\gcs_tables'  # Stores output GCS tables

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if wetted_dir == '':
        wetted_folder = dem_dir + '\\wetted_polygons\\small_increments'

    spatial_ref = arcpy.Describe(detrended_dem).spatialReference
    unit = spatial_ref.linearUnitName
    if unit == 'Meter':
        u = 'm'

    else:
        u = 'ft'

    for z in zs:
        z_str = float_keyz_format(z)
        label = z_str + u
        in_list = [wetted_dir + '\\wetted_poly_%sft.shp' % label, temp_files + '\\%scenterline_XS.shp' % label,
                   lines_dir + '\\%scenterline.shp' % label]

        if clip_poly != '' and os.path.exists(
                clip_poly):  # Allows a new/updated clip file to clip all data inputs and outputs, and create new XS for the clipped centerlines
            for j, file in enumerate(in_list):
                no_clip_name = file[:-4] + '_delete.shp'
                if os.path.exists(no_clip_name):
                    file_functions.delete_gis_files(no_clip_name)
                try:
                    arcpy.Rename_management(file, no_clip_name)
                    del_files.append(no_clip_name)
                except:
                    print('Error occured, could not rename %s file likely because it does not exist or is open' % file)

                if j != 1:
                    arcpy.Clip_analysis(no_clip_name, clip_poly, out_feature_class=file)

            create_station_lines_function(in_list[2], spacing, xs_lengths[loc_stage_index], stage=loc_stage)

        clipped_xs_lines = temp_files + '\\clipped_XS_lines_%s.shp' % label
        arcpy.Clip_analysis(in_list[1], in_list[0], out_feature_class=clipped_xs_lines)

        width_poly_loc = lines_dir + '\\width_rectangles_%s.shp' % label
        arcpy.Buffer_analysis(clipped_xs_lines, width_poly_loc, float(spacing / 2), line_side='FULL',
                              line_end_type='FLAT')
        arcpy.AddField_management(width_poly_loc, "Width", field_type="FLOAT")
        expression = ("(float(!Shape.area!)) / %d" % spacing)
        arcpy.CalculateField_management(width_poly_loc, "Width", expression, "PYTHON3")
        print('Width polygons for %sft stage created...' % z)

        arcpy.AddField_management(width_poly_loc, field_name="loc_id", field_type="SHORT")
        field_calc = "(int(!LOCATION!))"
        arcpy.CalculateField_management(width_poly_loc, field="loc_id", expression=field_calc,
                                        expression_type="PYTHON3")
        zonal_table = arcpy.sa.ZonalStatisticsAsTable(width_poly_loc, "loc_id", detrended_dem,
                                                      out_table=(temp_files + '\\stats_table_%s.dbf' % label),
                                                      statistics_type="ALL")
        width_poly = arcpy.JoinField_management(width_poly_loc, "loc_id", join_table=zonal_table, join_field="loc_id",
                                                fields=["MEAN"])

        csv_loc = out_dir + "\\%s_gcs_table.csv" % label
        tableToCSV(width_poly, csv_filepath=csv_loc, fld_to_remove_override=[])
        df = pd.read_csv(csv_loc)
        df.rename({'LOCATION': 'dist_down', 'Width': 'W', 'MEAN': 'Z'}, axis=1, inplace=True)
        df.sort_values(by=['dist_down'], inplace=True)
        df.to_csv(csv_loc)

        main_classify_landforms(tables=[csv_loc], w_field='W', z_field='Z', dist_field='dist_down', out_folder=out_dir)
        gcs_df = pd.read_csv(csv_loc)
        gcs_df.sort_values(by=['dist_down'], inplace=True)
        gcs_df.to_csv(csv_loc)
        print('GCS csv file made for stage %s...' % label)

    # Delete extra files if necessary
    for file in del_files:
        file_functions.delete_gis_files(file)

    print('GCS tables completed @ %s' % out_dir)
    return out_dir
