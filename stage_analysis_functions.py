def stage_analysis(detrended_dem, zs):
    """Runs if Stage Analysis is selected on the GCS Analysis window"""

    # make descriptive stats Excel file
    descriptive_stats_xlxs(detrended_dem, zs)
    return
