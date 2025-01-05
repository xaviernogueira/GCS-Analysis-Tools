import os
import logging

from processing_funcs import (
    prep_small_inc,
    pdf_cdf_plotting,
    stage_centerlines,
)

logger = logging.getLogger()

class WetterController:
    def __init__(self) -> None:
        self.maxs = [0]
        self.runs = 0
        self.highest = max(self.maxs)
        self.switch = False


def model_each_flow_stage(
    detrended_dem,
    max_stage,
    controller: WetterController,
) -> None:

    controller.runs += 1
    if max_stage > controller.highest:
        controller.switch = True

    out_dir = os.path.dirname(detrended_dem)
    out_list = []
    wet_dir = out_dir + '\\wetted_polygons'

    # if you run for the first time it starts from 0 to max stage, generating shapefiles and plotting
    # or the updated max_stage is higher, shapefiles are made for the missing stages, and plots are updated
    if controller.runs == 0 or controller.switch:
        logger.info('Making wetted area polygons...')
        prep_small_inc(
            detrended_dem=detrended_dem,
            max_stage=max_stage,
        )
        logger.info('Done')

        logger.info('Creating flow stage analysis plots...')
        imgs = pdf_cdf_plotting(
            in_dir=wet_dir,
            out_folder=out_dir,
            max_stage=max_stage,
        )
        logger.info('Done')
        logger.info('Plots @ %s' % os.path.dirname(imgs[0]))

        for img in imgs:
            name = os.path.basename(img)[:-4]
            out_list.append((name, img))

        controller.switch = False

    # if the max stage is less than on equal to the previous max stage, plots are updated
    else:
        logger.info('Creating flow stage analysis plots...')
        imgs = pdf_cdf_plotting(
            in_dir=wet_dir,
            out_folder=out_dir,
            max_stage=max_stage,
        )
        logger.info('Done')
        logger.info('Plots @ %s' % os.path.dirname(imgs[0]))

        for img in imgs:
            out_list.append(('Flow stage vs wetted area', img))

    # update list of used max stage values
    controller.maxs.append(max_stage)


