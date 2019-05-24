import configparser
import os
import logging

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def get_configuration(config_path=None):
    """
    Read .ini config file from given path
    """
    ref_path_config = configparser.ConfigParser()
    if config_path is None:
        log.info('config path not provided, use default config')
        ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    else:
        ref_path_config.read(config_path)
    return ref_path_config
