import sys
import logging
import os
from typing import List, Union
from pathlib import Path

from .stats_funcs import (
    sankey_chi_squared,
    violin_ttest,
)

from .plotting_funcs import (
    gcs_plotter,
    nested_landform_sankey,
)
from .nested_analysis_prep_funcs import prep_for_nesting_analysis

sys.path.append(str(Path(__file__).parent.parent))
from utils import prep_key_zs  

def nesting_analysis(
    detrended_dem: str,
    analysis_dir: str,
    zs: Union[str, List[Union[float, int]]],
) -> None:
    """Runs if Nesting Analysis is selected on the GCS Analysis window"""

    logging.info(f'Running GCS nested analysis on flow stages: {zs}')

    # get z labels for the key zs
    zs = prep_key_zs(zs)

    gcs_dir = os.path.dirname(detrended_dem) + '/gcs_tables'
    logging.info(
        f'Creating an aligned_gcs.csv table in {gcs_dir}, '
        f'which stores all flow stages: {zs}'
    )
    aligned_gcs_csv = prep_for_nesting_analysis(
       detrended_dem,
       zs=zs,
    )

    # NOTE: THIS IS JUST FOR TESTING SPEED
    aligned_gcs_csv = gcs_dir + '/aligned_gcs.csv'
    # TODO: REMOVE AFTER!

    logging.info(f'Done. Saved @ {aligned_gcs_csv}')

    logging.info('Making nested GCS line plots...')

    # TODO: make this work!
    nested_dir = gcs_plotter(
        detrended_dem,
        analysis_dir,
        zs,
        fields=['Ws', 'Zs', 'Ws_Zs'],
        aligned_csv=aligned_gcs_csv,
        together=True,
    )
    logging.info(f'Done. Plots saved @ {nested_dir}.')

    logging.info('Running landform transition Chi-Squared test...')
    sankey_csv = sankey_chi_squared(
        zs,
        aligned_gcs_csv,
        analysis_dir,
        detrended_dem,
    )
    logging.info(f'Done. Results @ {sankey_csv}')

    for state in [False, True]:
        logging.info(
            f'Creating Sankey plot visualizing flow stage transitions w/ ignore_normal={state}...'
        )
        out_html = nested_landform_sankey(
            detrended_dem,
            analysis_dir,
            zs,
            ignore_normal=state,
        )
        logging.info(f'Done. HTML plot saved @ {out_html}')

    logging.info('Running T-test to verify significance of violin analysis...')
    #violin_df = pd.DataFrame
    # TODO: combine into t-test
    violin_csv = violin_ttest(
        detrended_dem,
        analysis_dir,
        zs,
        thresh=0.50,
    )
    logging.info(f'Done. Results .csv saved @ {violin_csv}')

    logging.info(f'Done with nesting analysis. See results @ {nested_dir}')
