import logging
import os
import stats_functions
import plotting_functions
from typing import List, Union

#gcs_dir = dem_dir + '\\gcs_tables'


def nesting_analysis(
    detrended_dem: str,
    analysis_dir: str,
    zs: Union[str, List[float, int]],
):
    """Runs if Nesting Analysis is selected on the GCS Analysis window"""
    logging.info(f'Running GCS nested analysis on flow stages: {zs}')
    # TODO: fix the nesting analysis prep functions
