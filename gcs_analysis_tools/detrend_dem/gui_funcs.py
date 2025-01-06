import os
import sys
import logging
from pathlib import Path

from .processing_funcs import fit_params_txt, make_residual_plot, linear_fit_plot, \
    diagnostic_quick_plot, detrend_that_raster, linear_fit, prep_xl_file

sys.path.append(str(Path(__file__).parent.parent))

from utils import string_to_list

logger = logging.getLogger(__name__)

def make_xyz_plot(xyz_csv):
    """Plots thalweg elevation profile is plotted but not saved"""
    # set up directory for plots
    out_dir = os.path.dirname(xyz_csv) + '\\detrending_plots'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make arrays. Save and display (pop up) the thalweg elevation plot
    out_tuple = prep_xl_file(
        xyz_csv,
        in_columns=[
            'LOCATION',
            'POINT_X',
            'POINT_Y',
            'Value',
        ],
    )

    plot = diagnostic_quick_plot(
        location_np=out_tuple[0],
        z_np=out_tuple[1],
        out_dir=out_dir,
    )
    logging.info('Created thalweg elevation profile plot!')
    return plot

def make_fit_plots(xyz_csv, breakpoints):
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
    out_tuple = prep_xl_file(
        xyz_csv,
        in_columns=[
            'LOCATION',
            'POINT_X',
            'POINT_Y',
            'Value',
        ],
    )

    fit_out = linear_fit(
        location_np=out_tuple[0],
        z_np=out_tuple[1],
        xyz_table_loc=xyz_csv,
        bp_list=breakpoint_list,
    )

    # save and display (pop up) linear fit and residual plots. Generate txt
    fit_plot = linear_fit_plot(
        location_np=out_tuple[0],
        z_np=out_tuple[1],
        fit_params=fit_out[0],
        fit_np=fit_out[1],
        out_dir=out_dir,
    )

    res_plot = make_residual_plot(
        location_np=out_tuple[0],
        residual_np=fit_out[2],
        r2=fit_out[3],
        out_dir=out_dir,
    )

    txt = fit_params_txt(
        fit_params=fit_out[0],
        bp_list=breakpoint_list,
        out_dir=out_dir,
    )
    logging.info(
        'Text file listing linear piecewise fit components @ %s' % txt)

    return fit_plot, res_plot, breakpoint_list

def detrend(
        xyz_csv,
        in_dem,
        aoi_shp
) -> None:
    """Detrends the input raster based on .csv stored x, y, z, z_fit values for thalweg points"""

    logging.info('Detrending DEM...')
    detrended_dem = detrend_that_raster(
        xyz_csv=xyz_csv,
        in_dem=in_dem,
        aoi_shp=aoi_shp,
    )
    logging.info('Done')
    logging.info('Detrended DEM @ %s' % detrended_dem)




