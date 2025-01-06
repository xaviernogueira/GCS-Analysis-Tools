import logging
from .stats_funcs import (
    descriptive_stats_xlxs,
)
from .plotting_funcs import (
    gcs_plotter,
    heat_plotter,
    landform_pie_charts,
)
from typing import List, Union


def run_stage_analysis(
    detrended_dem: str,
    analysis_dir: str,
    zs: Union[str, List[Union[float, int]]],
):
    """Runs if Stage Analysis is selected on the GCS Analysis window"""

    # make descriptive stats Excel file
    logging.info(
        'Writing descriptive stats and WW-Runs test output to .xlsx...')
    out_xlsx = descriptive_stats_xlxs(
        zs=zs,
        analysis_dir=analysis_dir,
        detrended_dem=detrended_dem,
    )
    logging.info(f'Done! Output @ {out_xlsx}')

    # make stage based plots
    logging.info('Making GCS plots for each stage...')
    stage_plots_dir = gcs_plotter(
        detrended_dem=detrended_dem,
        analysis_dir=analysis_dir,
        zs=zs,
        together=False,
    )
    logging.info(f'Done! Output @ {stage_plots_dir}')

    # make Ws-Zs heatplots
    logging.info('Making GCS plots for each stage...')
    stage_plots_dir = heat_plotter(
        detrended_dem=detrended_dem,
        analysis_dir=analysis_dir,
        zs=zs,
        together=False,
    )
    logging.info(f'Done! Output @ {stage_plots_dir}')

    # make landform pie charts
    logging.info('Making GCS landforms pie charts for each stage...')
    stage_plots_dir = landform_pie_charts(
        detrended_dem=detrended_dem,
        analysis_dir=analysis_dir,
        zs=zs,
        together=False,
    )
    logging.info(f'Done! Output @ {stage_plots_dir}')
    logging.info('Stage based analysis complete.')
