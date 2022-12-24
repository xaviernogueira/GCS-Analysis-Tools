import arcpy
from stats_functions import *
from plotting_functions import *


def prep_locations(
    detrend_folder: str,
    flip: bool = False
) -> str:
    """This function takes a reach and creates a new csv with aligned location identifiers using a Thiessen Polygon methodology.
    Returns aligned_locations.csv in the \\landform_analysis sub-folder. This csv can be used to align any data field for any key z or stage range.
    When flip==True (false is default) locations are flipped to correctlyly link to an ALREADY FLIPPED GCS stage table. If neither table is flipped, use
    the flipped table function in plotting functions file!"""

    logging.info(
        'Creating aligned_locations.csv with aligned centerline locations / dist_down...')
    detrended_raster = detrend_folder + '\\ras_detren.tif'
    landform_folder = detrend_folder + '\\landform_analysis'
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

    arcpy.env.overwriteOutput = True
    arcpy.env.extent = detrended_raster
    del_files = []

    # TODO: get these functions working
    centerline_nums = find_centerline_nums(detrend_folder)
    spacing = find_xs_spacing(detrend_folder)

    for num in centerline_nums:
        line_loc = ('%s\\stage_centerline_%sft_DS.shp' %
                    (centerline_folder, num))
        station_lines = create_station_lines_function(
            line_loc,
            spacing=spacing,
            xs_length=5,
            stage=[],
        )
        station_lines = centerline_folder + \
            ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (num, spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis(
            [station_lines, line_loc],
            out_feature_class=(centerline_folder +
                               "\\station_points_%sft.shp" % num),
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
                except KeyError:
                    logging.warning('Cant delete field: %s' % field)
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
                centerline_folder +
                "\\station_points_%sft.shp" % min(centerline_nums),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )
        elif counter > 1:
            arcpy.Identity_analysis(
                centerline_folder + ("\\align_points%s.shp" %
                                     (int(counter - 1))),
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )

      # Creates a csv with the aligned locations for each centerline. Allows joins to add any data to this for analysis
    index_field = 'loc_%sft' % min(centerline_nums)
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


def nesting_analysis():
    """Runs if Nesting Analysis is selected on the GCS Analysis window"""
    return
