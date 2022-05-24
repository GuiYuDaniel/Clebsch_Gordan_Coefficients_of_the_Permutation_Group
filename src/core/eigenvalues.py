# -*- coding: utf-8 -*-
"""
this code for creating eigenvalues of 2-cycle class by Sn, young_diagrams
and provide a function to calc eigenvalue2yds
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第5节：置换群的CSCO-Ⅱ
# 这里计算二循环类算符是为了标记寻找唯一的本征值集"航道"

# 公式：
# λ = 1/2 * sum{l}(ν_l * (ν_l - 2l + 1))
#   = 1/2 * sum{i}(ν_i * (ν_i - 2l - 1))
# 其中，l是杨盘图的行号（从1开始），ν_l是指第l行的格子数
# 因为python中，index从0开始，所以l = i + 1


import copy
import time
from conf.cgc_config import default_s_n
from core.characters_and_gi import load_characters_and_gi
from core.young_diagrams import load_young_diagrams, is_young_diagram
from core.cgc_utils.cgc_db_typing import CGOrderInfo
from core.cgc_utils.cgc_local_db import get_cg_order_file_name, get_cg_order_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


def _calc_single_eigenvalue_of_2_cycle_class_of_yd(yd: list):
    """
    公式：
    λ = 1/2 * sum{l}(ν_l * (ν_l - 2l + 1))
      = 1/2 * sum{i}(ν_i * (ν_i - 2l - 1))
    其中，l是杨盘图的行号（从1开始），ν_l是指第l行的格子数
    因为python中，index从0开始，所以l = i + 1
    """
    # TODO check with source code
    if not is_young_diagram(yd):
        err_msg = "cannot calc eigenvalue with wrong young_diagram={}".format(yd)
        logger.error(err_msg)
        return False, err_msg
    if yd == [1]:
        return True, 1  # 注意，二循环类算符，本应从S2开始，其S1的值不是计算的，而是定义的延拓

    rst = 0
    for row in range(len(yd)):
        # rst_part = yd[row] * (yd[row] + 1 - 2 * (row + 1)) / 2
        rst_part = yd[row] * (yd[row] - 2 * row - 1)  # /2放外面
        rst += rst_part
    rst = int(rst / 2)
    return True, rst


def calc_eigenvalues_of_2_cycle_class_of_s_n(yd_list: list):
    if not isinstance(yd_list, list):
        err_msg = "yd_list={} for calc_eigenvalues_of_2_cycle_class must be list but not".format(yd_list)
        logger.error(err_msg)
        return False, err_msg

    eigenvalues_list = []
    for single_yd in yd_list:
        flag, single_eigenvalue = _calc_single_eigenvalue_of_2_cycle_class_of_yd(single_yd)
        if not flag:
            err_msg = "_calc_single_eigenvalue_of_2_cycle_class_of_yd meet error with single_yd={}, msg={}".format(
                single_yd, single_eigenvalue)
            logger.error(err_msg)
            return False, err_msg
        eigenvalues_list.append(single_eigenvalue)

    return True, eigenvalues_list


def get_yd_list_by_eigenvalue(s_n: int, eigenvalue: int):
    pass

