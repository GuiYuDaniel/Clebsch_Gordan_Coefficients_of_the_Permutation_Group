# -*- coding:utf8 -*-


import os
from conf.cgc_config import cgc_rst_folder
from utils.config import singleton_config
from utils.log import get_logger


logger = get_logger(__name__)


def _create_results_folder():
    top_path = singleton_config.top_path
    sys_rst_name = singleton_config.result_folder
    cgc_rst_full_path = os.path.join(top_path, sys_rst_name, cgc_rst_folder)
    if not isinstance(cgc_rst_full_path, str):
        err_msg = "results folder name={} type must be str but not, pls check".format(cgc_rst_full_path)
        logger.error(err_msg)
        raise Exception(err_msg)
    if not os.path.exists(cgc_rst_full_path):
        logger.info("init cgc results folder as {}".format(cgc_rst_full_path))
        os.makedirs(cgc_rst_full_path)
    return cgc_rst_full_path
