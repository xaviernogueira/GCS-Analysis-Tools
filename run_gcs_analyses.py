from stats_functions import *
from plotting_functions import *


def stage_analysis(detrended_dem, zs):
    """Runs if Stage Analysis is selected on the GCS Analysis window"""

    # make descriptive stats Excel file
    descriptive_stats_xlxs(detrended_dem, zs)
    return


def nesting_analysis():
    """Runs if Nesting Analysis is selected on the GCS Analysis window"""
    return