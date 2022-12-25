import logging
import os
import pathlib
from typing import List
import arcpy
import file_functions
from create_station_lines import create_station_lines_function


def find_xs_length(
        dem_dir: str,
        zs: List[float]):
    """This function takes a detrend folder location for a given reacg, as well as a list containing centerline_nums, and by using
    string splicing and Arc geomoetry objects returns a list containing the XS widths for each centerline_num XS file"""

    centerline_folder = dem_dir + '\\analysis_centerline_and_XS'
    full_list = [f for f in os.path.listdir(centerline_folder) if (
        f[21:26] == 'DS_XS' or f[22:27] == 'DS_XS') and f[-4:] == '.shp']
    full_list = [i for i in full_list if i[-8:]
                 != 'x5ft.shp' and i[-10:] != 'delete.shp']
    xs_lengths = []

    for z in zs:
        if os.path.exists(centerline_folder + '\\stage_centerline_%sft_DS_XS_%sft.shp' % (num, find_xs_spacing(detrend_folder))):
            xs_file = centerline_folder + \
                '\\stage_centerline_%sft_DS_XS_%sft.shp' % (
                    num, find_xs_spacing(detrend_folder))
        else:
            print(
                'Cant find XS file, make sure formatting is: stage_centerline_[CENTERLINE_NUM]ft_DS_XS_[XS SPACING]ft.shp')

        temp_list = []
        for row in arcpy.da.SearchCursor(xs_file, ["SHAPE@LENGTH"]):
            temp_list.append(int(row[0]))
        xs_lengths.append(temp_list[0])

    return xs_lengths


def find_xs_spacing(detrend_folder):
    """This function takes the detrended folder and centerline_nums (optional)"""
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    xs_files = [i for i in os.path.listdir(centerline_folder) if
                os.path.isfile(os.path.join(centerline_folder, i)) and i[-5:] == 't.shp' and len(i) > 32]
    xs_files = [i for i in xs_files if i[-8:] != 'x5ft.shp']

    temp_location = []
    cursor = arcpy.SearchCursor(centerline_folder + '\\%s' % xs_files[0])
    for row in cursor:
        temp_location.append(int(row.getValue('LOCATION')))
    temp_location.sort()
    spacing1 = int(temp_location[1] - temp_location[0])

    if xs_files[0][-8] == '_':
        spacing2 = int(xs_files[0][-7])
    elif xs_files[0][-9] == '_':
        spacing2 = int(xs_files[0][-8:-6])

    if spacing1 == spacing2:
        print('XS spacing is %sft...' % spacing1)
        return spacing1

    else:
        print('XS shapefile name spacing (from string splicing) is not equal to spacing found via arcpy Search Cursor.')
        return spacing2


def prep_for_nesting_analysis(
    detrended_dem: str,
    flip: bool = False
) -> str:
    """This function takes a reach and creates a new csv with aligned location identifiers using a Thiessen Polygon methodology.
    Returns aligned_locations.csv in the \\landform_analysis sub-folder. This csv can be used to align any data field for any key z or stage range.
    When flip==True (false is default) locations are flipped to correctlyly link to an ALREADY FLIPPED GCS stage table. If neither table is flipped, use
    the flipped table function in plotting functions file!"""
    # import arcpy only if necessary (should only need to be ran once)

    logging.info(
        'Creating aligned_locations.csv with aligned centerline locations / dist_down...')

    # set up env variables
    arcpy.env.overwriteOutput = True
    arcpy.env.extent = detrended_dem
    del_files = []

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    lines_dir = dem_dir + '\\centerlines'
    wetted_dir = dem_dir + '\\wetted_polygons'
    temp_files = dem_dir + '\\temp_files'
    out_dir = dem_dir + '\\gcs_tables'  # Stores output GCS tables

    # TODO: get these functions working
    # TODO: zs is basically Zs, but we do need their spacing
    spacing = find_xs_spacing(dem_dir)

    for z in zs:
        line_loc = ('%s\\stage_centerline_%sft_DS.shp' %
                    (lines_dir, z))
        station_lines = create_station_lines_function(
            line_loc,
            spacing=spacing,
            xs_length=5,
            stage=[],
        )
        station_lines = lines_dir + \
            ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (z, spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis(
            [station_lines, line_loc],
            out_feature_class=(lines_dir +
                               "\\station_points_%sft.shp" % z),
            join_attributes="ALL",
            output_type="POINT",
        )

        if z != min(zs):
            theis_loc = lines_dir + "\\thiessen_%sft.shp" % z
            arcpy.CreateThiessenPolygons_analysis(
                station_points,
                theis_loc,
                'ALL',
            )
            arcpy.AddField_management(
                theis_loc,
                ('loc_%sft' % z),
                'SHORT',
            )
            arcpy.CalculateField_management(
                theis_loc,
                ('loc_%sft' % z),
                expression='!LOCATION!',
                expression_type='PYTHON3',
            )

            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [('loc_%sft' % z), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except KeyError:
                    logging.warning('Cant delete field: %s' % field)
            arcpy.DeleteField_management(
                theis_loc,
                del_fields,
            )

    max_count = 0
    for counter, num in enumerate(zs):
        theis_loc = (lines_dir + "\\thiessen_%sft.shp" % num)
        out_points = lines_dir + ("\\align_points%s.shp" % counter)
        del_files.append(theis_loc)
        del_files.append(out_points)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(
                lines_dir +
                "\\station_points_%sft.shp" % min(zs),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )
        elif counter > 1:
            arcpy.Identity_analysis(
                lines_dir + ("\\align_points%s.shp" %
                             (int(counter - 1))),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )

      # Creates a csv with the aligned locations for each centerline. Allows joins to add any data to this for analysis
    index_field = 'loc_%sft' % min(zs)
    aligned_csv = landform_folder + '\\aligned_locations.csv'

    aligned_df = pd.read_csv(
        file_functions.table_to_csv(
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
        loc_fields = [j for j in list(
            out_aligned_df.columns.values) if j[:3] == 'loc']
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

    logging.info('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    logging.info('Empty aligned csv created @ %s!' % aligned_csv)
    return aligned_csv
