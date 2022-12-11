import arcpy
from arcpy import da
import file_functions
from file_functions import *
import create_station_lines
import statistics
import pandas as pd
import os
from typing import List, Union

arcpy.env.overwriteOutput = True


@err_info
def clean_in_table(
    table: pd.DataFrame,
    w_field: str = 'W',
    z_field: str = 'Z',
    dist_field: str = 'dist_down',
) -> pd.DataFrame:
    """Renames DataFrame columns corresponding to W, Z, and dist_down, if they are not already columns"""
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

def standardize(
    table: pd.DataFrame,
    fields: List[str]
) -> pd.DataFrame:
    """Makes standardized version of field in csv table by subtracting each value by mean and dividing by standard deviation.
    Creates a new column Ws_Zs storing C(Ws,Zs) values"""
    
    check_use(table)
    df = pd.read_csv(table)

    s_fields = []
    for f in fields:
        new_field = f + 's'
        df[new_field] = (df[f] - np.mean(df[f])) * 1.0 / np.std(df[f])
        s_fields.append(new_field)

    df['%s_%s' % (s_fields[0], s_fields[1])] = df[s_fields[0]] * df[s_fields[1]]
    
    return df.to_csv(
        table,
        index=False,
    )

def landforms(
    table: pd.DataFrame,
    zs_field: str = 'Zs',
    ws_field: str = 'Ws',
    na: int = -9999
) -> pd.DataFrame:
    """Classifies each row by corresponding landform type:
        oversized, nozzle, constricted pool, wide bar, normal channel
        Adds columns to input table"""

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


def main_classify_landforms(
    tables: List[str],
    w_field: str,
    z_field: str,
    dist_field: str
):
    """Classifies river segments as normal, wide bar, constricted pool, oversized, or nozzle

    Args:
        tables: a list of attribute table filenames for each set of wetted polygon rectangular XS's
        w_field: name of the attribute table field corresponding to width
        z_field: name of the attribute table field corresponding to detrended bed elevation
        dist_field: name of the attribute table field corresponding to distance downstream

    Returns:
        For each input table:
            a .csv containing dist_down, W, Z, Ws, Zs, Ws_Zs, and landform classification/code fields
            adds these computed values to attribute tables of corresponding wetted polygon rectangular XS's
    """
    logging.info('Classifying landforms...')
    #TODO: verify this behavior
    out_polys = []
    fields = [w_field, z_field]

    for i in range(len(tables)):
        table = tables[i]
        clean_in_table(table, w_field=w_field, z_field=z_field, dist_field=dist_field)
        standardize(table, fields=fields)
        landforms(table)

    logging.info('Finished.')

    return out_polys


# Main function that conducts GCS analysis
############################################


def extract_gcs(
    detrended_dem: str,
    zs: Union[str, List[Union[int, float]]],
    xs_lengths: Union[str, List[Union[int, float]]],
    spacing: Union[str, List[Union[int, float]]],
    clip_poly: str = ''
) -> List[str]:
    '''This function does a full GCS analysis using three specific key Zs that can include any float. Results saved
    to the gcs_ready_tables, as well as plotted. Results are aligned to the existing csv to facilitate landform analysis
    detrend
    wetted_folder is the folder containing small increment wetted polygons
    If key_zs parameter is an empty list, a range from 0 to max_stage (deafult is 20) makes gcs csvs at 1ft increments.
    wetted_folder parameter (optional) allows for a specific folder containing wetted polygons to be used instead of the assumed file structure.'''

    arcpy.env.overwriteOutput = True

    del_files = []  # stored deletable files
    out_csvs = []  # stores output gcs csv files

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

    lines_dir = dem_dir + '\\centerlines'
    wetted_dir = dem_dir + '\\wetted_polygons'
    temp_files = dem_dir + '\\temp_files'
    out_dir = dem_dir + '\\gcs_tables'  # Stores output GCS tables

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    og_dem = dem_dir + '\\las_dem.tif'

    # Get units for string labeling
    u, unit, spatial_ref = file_functions.get_label_units(detrended_dem)

    for i, z in enumerate(zs):
        z_str = float_keyz_format(z)
        label = z_str + u
        in_list = [
            wetted_dir + '\\wetted_poly_%s.shp' % label,
            temp_files + '\\%s_centerline_XS.shp' % label,
            lines_dir + '\\%s_centerline.shp' % label,
            ]

        xs_length = xs_lengths[i]

        # Allows a new/updated clip file to clip all data inputs and outputs, and create new XS for the clipped centerlines
        if clip_poly != '' and os.path.exists(clip_poly):
            for j, file in enumerate(in_list):
                no_clip_name = file.replace('.shp', '_delete.shp')
                if os.path.exists(no_clip_name):
                    arcpy.Delete_management(no_clip_name)
                try:
                    arcpy.Rename_management(file, no_clip_name)
                    del_files.append(no_clip_name)
                except PermissionError:
                    print('Permission Error: Could not rename %s file likely because it does not exist or is open' % file)

                if j != 1:
                    arcpy.Clip_analysis(
                        no_clip_name,
                        clip_poly,
                        out_feature_class=file,
                    )

        create_station_lines.create_station_lines_function(
            in_list[2],
            spacing,
            xs_length,
        )

        # Clip cross-sections by wetted area and create width rectangles
        clipped_xs_lines = temp_files + '\\clipped_XS_lines_%s.shp' % label
        arcpy.Clip_analysis(
            in_list[1],
            in_list[0],
            out_feature_class=clipped_xs_lines,
        )

        width_poly_loc = lines_dir + '\\width_rectangles_%s.shp' % label
        arcpy.Buffer_analysis(
            clipped_xs_lines,
            width_poly_loc,
            float(spacing / 2),
            line_side='FULL',
            line_end_type='FLAT',
        )

        # Extract rectangle width values and create a new attribute table Width field
        arcpy.AddField_management(
            width_poly_loc,
            "Width",
            field_type="FLOAT",
        )

        expression = ("(float(!Shape.area!)) / %d" % spacing)
        arcpy.CalculateField_management(
            width_poly_loc,
            "Width",
            expression,
            "PYTHON3",
        )
        print('Width polygons for %sft stage created...' % z)

        arcpy.AddField_management(
            width_poly_loc,
            field_name="loc_id",
            field_type="SHORT",
        )

        field_calc = "(int(!LOCATION!))"
        arcpy.CalculateField_management(
            width_poly_loc,
            field="loc_id",
            expression=field_calc,
            expression_type="PYTHON3",
        )

        # Extract the mean detrended DEM raster value from below each width rectangle and join to the shapefile
        zonal_table = arcpy.sa.ZonalStatisticsAsTable(
            width_poly_loc,
            "loc_id",
            detrended_dem,
            out_table=(temp_files + '\\stats_table_%s.dbf' % label),
            statistics_type="ALL",
        )

        arcpy.JoinField_management(
            width_poly_loc,
            "loc_id",
            join_table=zonal_table,
            join_field="loc_id",
            fields=["MEAN"],
        )

        # If las_dem.tif is unmoved and not renamed, we can extract the MEDIAN value and use it to flip tables if necessary
        zonal_table = arcpy.sa.ZonalStatisticsAsTable(
            width_poly_loc,
            "loc_id",
            og_dem,
            out_table=(temp_files + '\\no_detrend_stats_table_%s.dbf' % label),
            statistics_type="ALL",
        )

        no_sort = arcpy.JoinField_management(
            width_poly_loc,
            "loc_id",
            join_table=zonal_table,
            join_field="loc_id",
            fields=["MIN"],
        )

        # Sort width polygon by location identifier
        arcpy.Sort_management(
            no_sort,
            width_poly_loc.replace('.shp', '_s.shp'),
            [["LOCATION", "Ascending"]],
        )
        arcpy.Delete_management(no_sort)
        
        width_poly = arcpy.Rename_management(
            width_poly_loc.replace('.shp', '_s.shp'),
            width_poly_loc,
        )

        # Create flipped dist_down index if necessary
        cursor = arcpy.SearchCursor(width_poly)
        elevs = []
        locations = []
        for row in cursor:
            elevs.append(row.getValue("MIN"))
            locations.append(row.getValue("LOCATION"))
        arcpy.AddField_management(
            width_poly,
            'dist_down',
            "LONG",
        )
        max_loc = int(max(locations))

        try:
            if statistics.mean(elevs[10:20]) < statistics.mean(elevs[-20:-10]):
                expression = ("%d - (int(!LOCATION!))" % max_loc)
                arcpy.CalculateField_management(
                    width_poly,
                    'dist_down',
                    expression,
                    "PYTHON3",
                )

            else:
                expression = "(int(!LOCATION!))"
                arcpy.CalculateField_management(
                    width_poly,
                    'dist_down',
                    expression,
                    "PYTHON3",
                )
        except IndexError: raise ValueError(
            'Error: Cross-section series not longer than 20, cant establish which side is upstream.'
            )

        # Convert width polygon attribute table to a csv and classify landforms
        csv_loc = out_dir + "\\%s_gcs_table.csv" % label
        tableToCSV(
            width_poly,
            csv_filepath=csv_loc,
            fld_to_remove_override=[],
        )
        df = pd.read_csv(csv_loc)
        df.rename(
            {
                'LOCATION': 'location',
                'Width': 'W',
                'MEAN': 'Z',
                'MIN': 'Z_no_detrend',
            },
            axis=1,
            inplace=True,
        )
        df.sort_values(
            by=['dist_down'],
            inplace=True,
        )
        df = df[['location', 'dist_down', 'Z_no_detrend', 'W', 'Z']]
        df.to_csv(csv_loc)

        # calculate Ws, Zs, Ws_Zs, and landform codes
        main_classify_landforms(
            tables=[csv_loc],
            w_field='W',
            z_field='Z',
            dist_field='dist_down',
        )
        gcs_df = pd.read_csv(csv_loc)
        gcs_df.to_csv(csv_loc)
        print('GCS csv file made for stage %s...' % label)
        out_csvs.append(csv_loc)

    for file in del_files:
        file_functions.delete_gis_files(file)

    print('GCS tables completed @ %s' % out_dir)
    return out_csvs


def prep_locations(
    detrend_folder: str,
    flip: bool = False
) -> str:
    '''This function takes a reach and creates a new csv with aligned location identifiers using a Thiessen Polygon methodology.
    Returns aligned_locations.csv in the \\landform_analysis sub-folder. This csv can be used to align any data field for any key z or stage range.
    When flip==True (false is default) locations are flipped to correctlyly link to an ALREADY FLIPPED GCS stage table. If neither table is flipped, use
    the flipped table function in plotting functions file!'''

    print('Creating aligned_locations.csv with aligned centerline locations / dist_down...')
    detrended_raster = detrend_folder + '\\ras_detren.tif'
    landform_folder = detrend_folder + '\\landform_analysis'
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

    arcpy.env.overwriteOutput = True
    arcpy.env.extent = detrended_raster
    del_files = []

    #TODO: get these functions working
    centerline_nums = find_centerline_nums(detrend_folder)
    spacing = find_xs_spacing(detrend_folder)

    for num in centerline_nums:
        line_loc = ('%s\\stage_centerline_%sft_DS.shp' % (centerline_folder, num))
        station_lines = create_station_lines.create_station_lines_function(
            line_loc,
            spacing=spacing,
            xs_length=5,
            stage=[],
        )
        station_lines = centerline_folder + ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (num, spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis(
            [station_lines, line_loc],
            out_feature_class=(centerline_folder + "\\station_points_%sft.shp" % num),
            join_attributes="ALL",
            output_type="POINT",
        )

        if num != min(centerline_nums):
            theis_loc = centerline_folder + "\\thiessen_%sft.shp" % num
            arcpy.CreateThiessenPolygons_analysis(
                station_points,
                theis_loc,
                'ALL',
            )
            arcpy.AddField_management(
                theis_loc,
                ('loc_%sft' % num),
                'SHORT',
            )
            arcpy.CalculateField_management(
                theis_loc,
                ('loc_%sft' % num),
                expression='!LOCATION!',
                expression_type='PYTHON3',
            )

            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [('loc_%sft' % num), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except KeyError: print("Can't delete field: %s" % field)
            arcpy.DeleteField_management(
                theis_loc,
                del_fields,
            )

    max_count = 0
    for counter, num in enumerate(centerline_nums):
        theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
        out_points = centerline_folder + ("\\align_points%s.shp" % counter)
        del_files.append(theis_loc)
        del_files.append(out_points)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(
                centerline_folder + "\\station_points_%sft.shp" % min(centerline_nums),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )
        elif counter > 1:
            arcpy.Identity_analysis(
                centerline_folder + ("\\align_points%s.shp" % (int(counter - 1))),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )

      # Creates a csv with the aligned locations for each centerline. Allows joins to add any data to this for analysis
    index_field = 'loc_%sft' % min(centerline_nums)
    aligned_csv = landform_folder + '\\aligned_locations.csv'

    aligned_df = pd.read_csv(
        file_functions.tableToCSV(
            out_points,
            csv_filepath=aligned_csv,
            fld_to_remove_override=['FID_statio', 'FID_thiess'],
        ),
    )

    aligned_df.rename(
        columns={'LOCATION': index_field},
        inplace=True,
    )
    aligned_df.drop_duplicates(subset=[index_field], inplace=True)

    headers = list(aligned_df.columns.values)
    keep_headers = [i for i in headers if i[:3] == 'loc']

    out_aligned_df = aligned_df.loc[:, keep_headers]

    out_aligned_df.sort_values(
        by=[index_field],
        inplace=True,
    )
    out_aligned_df.set_index(
        index_field,
        inplace=True,
    )

    if flip:
        loc_fields = [j for j in list(out_aligned_df.columns.values) if j[:3] == 'loc']
        loc_nums = []

        for loc_field in loc_fields:
            if loc_field[5] == 'f':
                loc_nums.append(loc_field[4])
            else:
                loc_nums.append(loc_field[4:6])
            
            temp_max = np.nanmax(out_aligned_df.loc[:, loc_field].to_numpy())
            dist_list = out_aligned_df.loc[:, [loc_field]].squeeze().to_list()

            loc_np = np.array([int(temp_max - i) for i in dist_list])
            out_aligned_df[loc_field] = loc_np

        min_loc = loc_fields[loc_nums.index(min(loc_nums, key=int))]
        out_aligned_df.sort_values(
            str(min_loc),
            inplace=True,
        )
        out_aligned_df.to_csv(aligned_csv)

    else:
        out_aligned_df.to_csv(aligned_csv)

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    print('Empty aligned csv created @ %s!' % aligned_csv)
    return aligned_csv
