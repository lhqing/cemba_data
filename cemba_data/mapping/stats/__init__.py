import pathlib

from .mc import snmc_mapping_stats
from ...utilities import get_configuration


def mapping_stats(output_dir):
    """This is UID level mapping summary, the config file is in parent dir"""
    output_dir = pathlib.Path(output_dir).absolute()
    config = get_configuration(output_dir.parent / 'mapping_config.ini')
    mode = config['mode']

    if mode == 'mc':
        final_df = snmc_mapping_stats(output_dir, config)
    else:
        raise ValueError

    # save
    final_df.to_csv(output_dir / 'MappingSummary.csv.gz')
    return