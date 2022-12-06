# -*- coding:utf8 -*-


import os
from conf.cgc_config import cgc_db_name
from utils.config import singleton_config
from utils.log import get_logger


logger = get_logger(__name__)


def _create_results_folder():
    result_folder = singleton_config.result_folder
    if not os.path.exists(result_folder):
        raise FileExistsError("can not find result_folder={}, pls check".format(result_folder))
    cgc_rst_full_path = os.path.join(result_folder, cgc_db_name)
    if not os.path.exists(cgc_rst_full_path):
        logger.info("init cgc results folder as {}".format(cgc_rst_full_path))
        os.makedirs(cgc_rst_full_path)
    return cgc_rst_full_path


def concat_s_n(*args, s_n_type=int, **kwargs):
    """
    比对所有输入的s_n。如果一致，输出s_n
    """
    if not args and not kwargs:
        err_msg = "both args={} and kwargs={} must have one real".format(args, kwargs)
        logger.error(err_msg)
        return False, err_msg

    all_params_set = set(list(args) + list(kwargs.values()))
    if len(all_params_set) != 1:
        err_msg = "len of set(all_params_list) must eq 1 but not, with set(all_params_list)={} args={}, kwargs={}" \
                  "".format(all_params_set, args, kwargs)
        logger.error(err_msg)
        return False, err_msg

    rst = list(all_params_set)[0]
    if not isinstance(rst, s_n_type):
        err_msg = "concat rst={} should eq with s_n_type={} but not".format(rst, s_n_type)
        logger.error(err_msg)
        return False, err_msg

    return True, rst
