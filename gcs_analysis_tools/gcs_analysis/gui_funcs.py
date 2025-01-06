import os
import logging
from typing import Union, List

from .calculate_gcs_functs import extract_gcs
from .stage_analysis_funcs import run_stage_analysis
from .nesting_analysis_funcs import nesting_analysis


def run_gcs_analyses(
    detrended_dem: str,
    zs: Union[str, List[Union[int, float]]],
    xs_lengths: Union[str, List[Union[int, float]]],
    xs_spacing: Union[str, List[int]],
    clip_poly: str = '',
    stage_plots: bool = False,
    nest_plots: bool = False,
    analysis_dir: str = None,
) -> None:
    """Controls GCS calculation and analysis behavior"""
    # get path to base outputs on

    if not stage_plots and not nest_plots:
        logging.info('Extract GCS series...')
        extract_gcs(
            detrended_dem,
            zs,
            xs_lengths,
            xs_spacing,
            clip_poly=clip_poly,
        )
        logging.info('Done')

    # control GCS analysis output locations
    else:
        default_dir = os.path.dirname(
            detrended_dem) + '\\GCS_analysis_outputs'
        if analysis_dir == '':
            analysis_dir = default_dir
            logging.info(f'Analysis outputs will be @ {default_dir}')
        elif not os.path.exists(analysis_dir):
            analysis_dir = default_dir
            logging.warning(
                f'Desired output location cannot be found! '
                f'Defaulting to {default_dir}'
            )
        else:
            logging.info(f'Analysis outputs will be @ {analysis_dir}')

    # run desired GCS analyses!
    if stage_plots and not nest_plots:
        logging.info('Running flow stage level GCS analyses...')
        run_stage_analysis(
            detrended_dem,
            analysis_dir,
            zs,
        )
        logging.info('Done')

    elif nest_plots and not stage_plots:
        logging.info('Running nested GCS analyses...')
        nesting_analysis(
            detrended_dem,
            analysis_dir,
            zs,
        )
        logging.info('Done')
    else:
        logging.info(
            'Running both flow stage and nesting GCS analyses...')
        run_stage_analysis(
            detrended_dem,
            analysis_dir,
            zs,
        )
        nesting_analysis(
            detrended_dem,
            analysis_dir,
            zs,
        )
        logging.info('Done')


