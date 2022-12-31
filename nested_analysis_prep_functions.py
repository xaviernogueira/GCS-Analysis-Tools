import logging
import os
import arcpy
import file_functions
import pandas as pd
import numpy as np
from scipy import stats
from typing import List, Dict, Union
from create_station_lines import create_station_lines_function


def find_xs_spacing(
    lines_dir: str,
    z_labels: List[str],
    unit_label: str,
) -> Dict[str, int]:
    """This function takes the detrended folder and centerline_nums (optional)"""

    # find width polygons for each flow stage
    width_polygons_dict = {}
    for z_str in z_labels:
        width_path = lines_dir + f'\\width_rectangles_{z_str}{unit_label}.shp'
        if os.path.exists(width_path):
            width_polygons_dict[z_str] = {
                'polygon_path': width_path,
            }
        else:
            raise KeyError(
                f'Cant find the following width polygon: {width_path}'
            )

    # use SearchCursor to see gaps between centerlines, take the most frequent
    for z_str, sub_dict in width_polygons_dict.items():
        poly_path = sub_dict['polygon_path']

        locs = []
        cursor = arcpy.SearchCursor(poly_path)
        for row in cursor:
            locs.append(float(row.getValue('LOCATION')))
        locs.sort()

        spacings = []
        i = 0
        while int(i + 1) < len(locs):
            spacings.append(int(locs[i + 1] - locs[i]))
            i += 1
        spacing = stats.mode(np.array(spacings)).mode[0]
        width_polygons_dict[z_str]['spacing'] = spacing
    return width_polygons_dict


def prep_for_nesting_analysis(
    detrended_dem: str,
    zs: Union[str, List[Union[float, int]]],
) -> str:
    """This function takes a reach and creates a new csv with aligned location identifiers using a Thiessen Polygon methodology.
    Returns aligned_locations.csv in the \\landform_analysis sub-folder. This csv can be used to align any data field for any key z or stage range.
    When flip==True (false is default) locations are flipped to correctlyly link to an ALREADY FLIPPED GCS stage table. If neither table is flipped, use
    the flipped table function in plotting functions file!"""

    logging.info(
        'Creating aligned_gcs.csv with aligned centerline locations / dist_down...')

    # pull in flow stages
    if isinstance(zs, str):
        zs = file_functions.string_to_list(zs, format='float')
    elif isinstance(zs, list):
        zs = [float(z) for z in zs]
    else:
        raise ValueError(
            'Key flow stage parameter input incorrectly. '
            'Please enter stage heights separated only by commas (i.e. 0.2,0.7,3.6)'
        )

    z_labels = [file_functions.float_keyz_format(z) for z in zs]

    # set up env variables
    arcpy.env.overwriteOutput = True
    arcpy.env.extent = detrended_dem
    del_files = []

    # set up directories
    dem_dir = os.path.dirname(detrended_dem)
    lines_dir = dem_dir + '\\centerlines'
    gcs_dir = dem_dir + '\\gcs_tables'
    temp_files = dem_dir + '\\temp_files'

    if not os.path.exists(temp_files):
        os.mkdir(temp_files)

    # find units
    u = file_functions.get_label_units(detrended_dem)[0]

    # find cross-section spacings
    spacings_dict = find_xs_spacing(
        lines_dir,
        z_labels,
        unit_label=u,
    )

    # make theissen polygons for each flow stage station points
    for i, z in enumerate(zs):
        z_str = z_labels[i]
        z_xs_spacing = int(spacings_dict[z_str]['spacing'])
        line_shp = lines_dir + f'//{z_str}{u}_centerline.shp'

        station_lines = create_station_lines_function(
            line_shp,
            spacing=z_xs_spacing,
            xs_length=5,
        )

        # keep track of temp output for deletion
        station_lines = temp_files + f'//{z_str}{u}_centerline_XS.shp'
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis(
            [station_lines, line_shp],
            out_feature_class=(
                temp_files + f"\\station_points_{z_str}{u}.shp"
            ),
            join_attributes="ALL",
            output_type="POINT",
        )

        if z != min(zs):
            theis_loc = temp_files + f"\\thiessen_{z_str}{u}.shp"
            arcpy.CreateThiessenPolygons_analysis(
                station_points,
                theis_loc,
                'ALL',
            )
            arcpy.AddField_management(
                theis_loc,
                (f'loc_{z_str}{u}'),
                'SHORT',
            )
            arcpy.CalculateField_management(
                theis_loc,
                (f'loc_{z_str}{u}'),
                expression='!LOCATION!',
                expression_type='PYTHON3',
            )

            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [(f'loc_{z_str}{u}'), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except KeyError or PermissionError:
                    logging.warning('Cant delete field: %s' % field)
            arcpy.DeleteField_management(
                theis_loc,
                del_fields,
            )

    # use identity analysis to find nearest baseline cross-section index for each stage
    max_count = 0
    min_z_str = file_functions.float_keyz_format(min(zs))
    for counter, z_str in enumerate(z_labels):
        theis_loc = temp_files + f"\\thiessen_{z_str}{u}.shp"
        out_points = temp_files + ("\\align_points%s.shp" % counter)
        del_files.append(theis_loc)
        del_files.append(out_points)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(
                temp_files +
                f"\\station_points_{min_z_str}{u}.shp",
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )
        elif counter > 1:
            c = int(counter - 1)
            arcpy.Identity_analysis(
                temp_files + f"\\align_points{c}.shp",
                theis_loc,
                out_feature_class=out_points,
                join_attributes='ALL',
            )

    # create aligned_gcs.csv storing all the data along the baseflow centerline index
    index_field = f'loc_{min_z_str}{u}'
    aligned_csv = gcs_dir + '\\aligned_gcs.csv'

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

    out_aligned_df = aligned_df[keep_headers]

    out_aligned_df.sort_values(
        by=[index_field],
        inplace=True,
    )

    # join the data for each Z to the dataframe
    for i, z in enumerate(zs):
        z_str = z_labels[i]
        join_col = f'loc_{z_str}{u}'
        stage_df = pd.read_csv(
            gcs_dir + f'\\{z_str}{u}_gcs_table.csv',
        ).set_index('dist_down')

        # delete excess columns
        del_cols = [str(i) for i in stage_df.columns if (
            str(i) == 'location') or ('Unnamed' in str(i))]

        stage_df = stage_df[[str(i)
                             for i in stage_df.columns if i not in del_cols]]

        # rename columns
        rename_dict = {}
        for col in stage_df.columns:
            col = str(col)
            rename_dict[col] = f'{z_str}{u}_{col}'
        stage_df.rename(columns=rename_dict, inplace=True)

        stage_df.index.name = join_col

        # join to out_aligned_df
        out_aligned_df.set_index(join_col, inplace=True)
        out_aligned_df = out_aligned_df.join(stage_df)
        out_aligned_df.reset_index(inplace=True)

    logging.info('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    # overwrite the previous csv
    out_aligned_df.to_csv(aligned_csv)
    logging.info('Aligned GCS .csv created @ %s!' % aligned_csv)

    return aligned_csv
