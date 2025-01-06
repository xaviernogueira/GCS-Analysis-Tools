import logging
import sys
from pathlib import Path

from .processing_funcs import river_builder_harmonics

sys.path.append(str(Path(__file__).parent.parent))
from utils import string_to_list

logger = logging.getLogger()

def export_to_river_builder(
    in_csv,
    index_field,
    units,
    r2_threshold,
    n_harmonics,
    methods,
    field_headers,
) -> None:
    """Cleans inputs and runs river builder export"""

    logger.info('Initiating export to River Builder...')

    # clean inputs
    in_csv = Path(in_csv)
    if not in_csv.exists():
        raise ValueError(
            f'Cant find input {str(in_csv)}'
        )
    else:
        in_csv = str(in_csv)

    if r2_threshold == '':
        r2_threshold = None
    else:
        r2_threshold = float(r2_threshold)

    if not r2_threshold > 0 and not r2_threshold <= 1:
        raise ValueError(
            f'Input R-Sqaured threshold is invalid! Must be > 0 but <= 1.'
        )

    if n_harmonics == '' or n_harmonics == '0':
        n_harmonics = None
    else:
        n_harmonics = int(n_harmonics)

    if field_headers == '':
        field_headers = None
    else:
        field_headers = string_to_list(field_headers, format='string')

    out_dir = river_builder_harmonics(
        in_csv,
        index_field,
        units,
        r2_threshold=r2_threshold,
        n_harmonics=n_harmonics,
        methods=methods,
        field_headers=field_headers,
    )
    logger.info(f'Done. Output @ {out_dir}')

