# -*- coding: utf-8 -*-
"""
this code for creating Characters and Gi of Permutation Groups by Sn
录入Sn的置换群特征标
拿来主义，不计算了!
"""

# 就是表示矩阵对角线上元素之和（trace迹）。对于有限群，双不表示的特征标正交归g
# 数据来源：
# 《The Theory of Group Characters and Matrix Representations of Groups》-2nd ed （Dudley E. Littlewood）1950 2006


import copy
import numpy as np
import time
from conf.cgc_config import default_s_n, min_s_n_of_characters_and_gi
from core.young_diagrams import load_young_diagrams
from core.cgc_utils.cgc_db_typing import CharacterAndGiInfo
from core.cgc_utils.cgc_local_db import get_characters_and_gi_file_name, get_characters_and_gi_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


class CharacterData(object):
    """Littlewood书中的特征标表"""

    def __init__(self):
        _ch_1 = {(1,): [1]}
        _ch_2 = {(2,): [1, 1],
                 (1, 1,): [1, -1]}
        _ch_3 = {(3,): [1, 1, 1],
                 (2, 1,): [2, 0, -1],
                 (1, 1, 1,): [1, -1, 1]}
        _ch_4 = {(4,): [1, 1, 1, 1, 1],
                 (3, 1,): [3, 1, 0, -1, -1],
                 (2, 2,): [2, 0, -1, 0, 2],
                 (2, 1, 1,): [3, -1, 0, 1, -1],
                 (1, 1, 1, 1,): [1, -1, 1, -1, 1]}
        _ch_5 = {(5,): [1, 1, 1, 1, 1, 1, 1],
                 (4, 1,): [4, 2, 1, 0, 0, -1, -1],
                 (3, 2,): [5, 1, -1, -1, 1, 1, 0],
                 (3, 1, 1,): [6, 0, 0, 0, -2, 0, 1],
                 (2, 2, 1,): [5, -1, -1, 1, 1, -1, 0],
                 (2, 1, 1, 1,): [4, -2, 1, 0, 0, 1, -1],
                 (1, 1, 1, 1, 1): [1, -1, 1, -1, 1, -1, 1]}
        _ch_6 = {(6,): [1] * 11,
                 (5, 1,): [5, 3, 2, 1, 1, 0, 0, -1, -1, -1, -1],
                 (4, 2,): [9, 3, 0, -1, 1, 0, -1, 0, 1, 3, 0],
                 (4, 1, 1,): [10, 2, 1, 0, -2, -1, 0, 1, 0, -2, 1],
                 (3, 3,): [5, 1, -1, -1, 1, 1, 0, 0, -1, -3, 2],
                 (3, 2, 1,): [16, 0, -2, 0, 0, 0, 1, 0, 0, 0, -2],
                 (2, 2, 2,): [5, -1, -1, 1, 1, -1, 0, 0, -1, 3, 2],
                 (3, 1, 1, 1,): [10, -2, 1, 0, -2, 1, 0, -1, 0, 2, 1],
                 (2, 2, 1, 1,): [9, -3, 0, 1, 1, 0, -1, 0, 1, -3, 0],
                 (2, 1, 1, 1, 1,): [5, -3, 2, -1, 1, 0, 0, 1, -1, 1, -1],
                 (1, 1, 1, 1, 1, 1,): [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1]}
        _ch_7 = {(7,): [1] * 15,
                 (6, 1,): [6, 4, 3, 2, 2, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1],
                 (5, 2,): [14, 6, 2, 0, 2, 0, -1, -1, 0, 2, -1, 1, 2, 0, 0],
                 (5, 1, 1,): [15, 5, 3, 1, -1, -1, 0, 0, -1, -3, 0, 0, -1, 1, 1],
                 (4, 3,): [14, 4, -1, -2, 2, 1, -1, 0, 0, 0, 2, -1, -1, 1, 0],
                 (4, 2, 1,): [35, 5, -1, -1, -1, -1, 0, 1, 1, 1, -1, 0, -1, -1, 0],
                 (3, 3, 1,): [21, 1, -3, -1, 1, 1, 1, 0, -1, -3, 0, 1, 1, -1, 0],
                 (4, 1, 1, 1,): [20, 0, 2, 0, -4, 0, 0, 0, 0, 0, 2, 0, 2, 0, -1],
                 (3, 2, 2,): [21, -1, -3, 1, 1, -1, 1, 0, -1, 3, 0, -1, 1, 1, 0],
                 (3, 2, 1, 1,): [35, -5, -1, 1, -1, 1, 0, -1, 1, -1, -1, 0, -1, 1, 0],
                 (2, 2, 2, 1,): [14, -4, -1, 2, 2, -1, -1, 0, 0, 0, 2, 1, -1, -1, 0],
                 (3, 1, 1, 1, 1,): [15, -5, 3, -1, -1, 1, 0, 0, -1, 3, 0, 0, -1, -1, 1],
                 (2, 2, 1, 1, 1,): [14, -6, 2, 0, 2, 0, -1, 1, 0, -2, -1, -1, 2, 0, 0],
                 (2, 1, 1, 1, 1, 1,): [6, -4, 3, -2, 2, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1],
                 (1, 1, 1, 1, 1, 1, 1,): [1, -1] * 7 + [1]}
        _ch_8 = {(8,): [1] * 22,
                 (7, 1,): [7, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0] + [-1] * 7,
                 (6, 2,): [20, 10, 5, 2, 4, 1, 0, -1, 0, 2, -1, 0, 1, -1, -1, 0, 0, 2, 1, 1, 0, 4],
                 (6, 1, 1,): [21, 9, 6, 3, 1, 0, 1, 0, -1, -3, 0, -1, -2, 0, 0, 1, 1, -1, 0, 0, 1, -3],
                 (5, 3,): [28, 10, 1, -2, 4, 1, -2, -1, 0, 2, 1, 0, 1, 1, 0, 0, 0, -2, -1, 1, 1, -4],
                 (5, 2, 1,): [64, 16, 4, 0, 0, -2, -1, 0, 0, 0, -2, 1, 0, 0, 1, 0, 0, 0, 0, -2, -1, 0],
                 (5, 1, 1, 1,): [35, 5, 5, 1, -5, -1, 0, 0, -1, -3, 2, 0, 1, 1, 0, -1, -1, 1, 0, 2, 0, 3],
                 (4, 4,): [14, 4, -1, -2, 2, 1, -1, 0, 0, 0, 2, -1, -1, 1, 0, 0, 2, 2, 0, -2, -1, 6],
                 (4, 3, 1,): [70, 10, -5, -4, 2, 1, 0, 1, 0, -2, 1, 0, -1, -1, 0, 0, -2, 0, 1, 1, 0, -2],
                 (4, 2, 2,): [56, 4, -4, 0, 0, -2, 1, 1, 0, 4, -1, -1] + [0] * 6 + [-1, 1, 1, 8],
                 (4, 2, 1, 1,): [90, 0, 0, 0, -6, 0, 0, 0, 2] + [0] * 5 + [-1, 0, 2, 0, 0, 0, 0, -6],
                 (3, 3, 2,): [42, 0, -6, 0, 2, 0, 2, 0, -2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, -1, -6],
                 (3, 3, 1, 1,): [56, -4, -4, 0, 0, 2, 1, -1, 0, -4, -1, 1] + [0] * 6 + [-1, -1, 1, 8],
                 (3, 2, 2, 1,): [70, -10, -5, 4, 2, -1, 0, -1, 0, 2, 1, 0, -1, 1, 0, 0, -2, 0, 1, -1, 0, -2],
                 (2, 2, 2, 2,): [14, -4, -1, 2, 2, -1, -1, 0, 0, 0, 2, 1, -1, -1, 0, 0, 2, -2, 0, 2, -1, 6],
                 (4, 1, 1, 1, 1,): [35, -5, 5, -1, -5, 1, 0, 0, -1, 3, 2, 0, 1, -1, 0, 1, -1, -1, 0, -2, 0, 3],
                 (3, 2, 1, 1, 1,): [64, -16, 4, 0, 0, 2, -1, 0, 0, 0, -2, -1, 0, 0, 1, 0, 0, 0, 0, 2, -1, 0],
                 (2, 2, 2, 1, 1,): [28, -10, 1, 2, 4, -1, -2, 1, 0, -2, 1, 0, 1, -1, 0, 0, 0, 2, -1, -1, 1, -4],
                 (3, 1, 1, 1, 1, 1,): [21, -9, 6, -3, 1, 0, 1, 0, -1, 3, 0, 1, -2, 0, 0, -1, 1, 1, 0, 0, 1, -3],
                 (2, 2, 1, 1, 1, 1,): [20, -10, 5, -2, 4, -1, 0, 1, 0, -2, -1, 0, 1, 1, -1, 0, 0, -2, 1, -1, 0, 4],
                 (2, 1, 1, 1, 1, 1, 1,): [7, -5, 4, -3, 3, -2, 2, -1, 1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1],
                 (1, 1, 1, 1, 1, 1, 1, 1,): [1, -1] * 10 + [1, 1]}
        # TODO 检查上面的(1, 1, 1, 1, 1, 1, 1, 1,) 最有有点不一样不是1，-1循环了
        _ch_9 = {(9,): [1] * 14 + [1] * 16,
                 (8, 1,): [8, 6, 5, 4, 4, 3, 3, 2, 2, 2, 2, 1, 1, 1] +
                          [1] + [0] * 7 + [-1] * 8,
                 (7, 2,): [27, 15, 9, 5, 7, 3, 2, 0, 1, 3, 0, 0, 1, -1,
                           -1, -1, -1, 1, 0, 0, -1, 3, 0, 0, 0, 1, 1, 2, 3, 0],
                 (7, 1, 1,): [28, 14, 10, 6, 4, 2, 3, 1, 0, -2, 1, -1, -2, 0,
                              0, 0, 0, -2, -1, -1, 0, -4, 1, 1, 1, 0, 0, -1, -2, 1],
                 (6, 3,): [48, 20, 6, 0, 8, 2, -2, -2, 0, 4, 0, 0, 2, 0,
                           -1, 0, 0, 0, 0, 2, 1, 0, 0, 3, 1, 0, -1, -2, -2, 0],
                 (6, 2, 1,): [105, 35, 15, 5, 5, -1, 0, -1, -1, -1, -3, 0, -1, -1,
                              0, 1, 1, 1, 1, -1, 0, 1, 0, -3, -1, -1, 0, 0, -1, 0],
                 (6, 1, 1, 1,): [56, 14, 11, 4, -4, -1, 1, 0, -2, -6, 2, -1, -1, 1,
                                 0, 0, 0, 0, 0, 2, 1, 0, -1, 2, 0, 1, 0, 1, 3, -1],
                 (5, 4,): [42, 14, 0, -4, 6, 2, -3, -1, 0, 2, 3, -1, 0, 2,
                           0, 0, 2, 0, -1, -1, 0, 2, 1, -3, -1, 0, 0, 1, 2, 0],
                 (5, 3, 1,): [162, 36, 0, -6, 6, 0, -3, 0, 0, 0, 0, 1, 0, 0,
                              1, 0, -2, -2, 0, 0, 0, -6, -1, 0, 0, 0, 1, 1, 0, 0],
                 (5, 2, 2,): [120, 20, 0, 0, 0, -4, 0, 1, 0, 4, -3, 0, 0, 0,
                              1, 0, 0, 0, -1, -1, 0, 8, 0, 3, 1, 0, -1, 0, 4, 0],
                 (5, 2, 1, 1,): [189, 21, 9, 1, -11, -3, -1, 0, 1, -3, 0, 1, 1, 1,
                                 0, -1, 1, 1, 0, 0, -1, -3, 1, 0, 0, 1, 0, -1, -3, 0],
                 (4, 4, 1,): [84, 14, -6, -6, 4, 2, -1, 1, 0, -2, 3, -1, -2, 0,
                              0, 0, 0, 2, 1, -1, -1, 4, -1, 3, 1, 0, 0, -1, -2, 0],
                 (4, 3, 2,): [168, 14, -15, -4, 4, -1, 3, 2, -2, 2, 0, -1, 1, -1,
                              0, 0, 0, 0, 0, 2, 0, 0, 1, -3, -1, 1, 0, -1, -1, 0],
                 (4, 3, 1, 1,): [216, 6, -9, -4, -4, 3, 1, 0, 2, -6, 0, 1, -1, -1,
                                 -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1, -1, 1, 3, 0],
                 (5, 1, 1, 1, 1,): [70, 0, 10, 0, -10, 0, 0, 0, -2, 0, 4, 0, 2, 0,
                                    0, 0, -2, 0, 0, 0, 0, 6, 0, -2, 0, -2, 0, 0, 0, 1],
                 (3, 3, 3,): [42, 0, -6, 0, 2, 0, 2, 0, -2, 0, 0, 0, 2, 0,
                              0, 0, 2, 0, 0, 0, -1, -6, 0, 6, 0, -2, 0, 2, 0, 0],
                 (4, 2, 2, 1,): [216, -6, -9, 4, -4, -3, 1, 0, 2, 6, 0, -1, -1, 1,
                                 -1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 1, 1, -3, 0],
                 (3, 3, 2, 1,): [168, -14, -15, 4, 4, 1, 3, -2, -2, -2, 0, 1, 1, 1,
                                 0, 0, 0, 0, 0, -2, 0, 0, -1, -3, 1, 1, 0, -1, 1, 0],
                 (3, 2, 2, 2,): [84, -14, -6, 6, 4, -2, -1, -1, 0, 2, 3, 1, -2, 0,
                                 0, 0, 0, -2, 1, 1, -1, 4, 1, 3, -1, 0, 0, -1, 2, 0],
                 (4, 2, 1, 1, 1,): [189, -21, 9, -1, -11, 3, -1, 0, 1, 3, 0, -1, 1, -1,
                                    0, 1, 1, -1, 0, 0, -1, -3, -1, 0, 0, 1, 0, -1, 3, 0],
                 (3, 3, 1, 1, 1,): [120, -20, 0, 0, 0, 4, 0, -1, 0, -4, -3, 0, 0, 0,
                                    1, 0, 0, 0, -1, 1, 0, 8, 0, 3, -1, 0, 1, 0, -4, 0],
                 (3, 2, 2, 1, 1,): [162, -36, 0, 6, 6, 0, -3, 0, 0, 0, 0, -1, 0, 0,
                                    1, 0, -2, 2, 0, 0, 0, -6, 1, 0, 0, 0, -1, 1, 0, 0],
                 (2, 2, 2, 2, 1,): [42, -14, 0, 4, 6, -2, -3, 1, 0, -2, 3, 1, 0, -2,
                                    0, 0, 2, 0, -1, 1, 0, 2, -1, -3, 1, 0, 0, 1, -2, 0],
                 (4, 1, 1, 1, 1, 1,): [56, -14, 11, -4, -4, 1, 1, 0, -2, 6, 2, 1, -1, -1,
                                       0, 0, 0, 0, 0, -2, 1, 0, 1, 2, 0, 1, 0, 1, -3, -1],
                 (3, 2, 1, 1, 1, 1,): [105, -35, 15, -5, 5, 1, 0, 1, -1, 1, -3, 0, -1, 1,
                                       0, -1, 1, -1, 1, 1, 0, 1, 0, -3, 1, -1, 0, 0, 1, 0],
                 (2, 2, 2, 1, 1, 1,): [48, -20, 6, 0, 8, -2, -2, 2, 0, -4, 0, 0, 2, 0,
                                       -1, 0, 0, 0, 0, -2, 1, 0, 0, 3, -1, 0, 1, -2, 2, 0],
                 (3, 1, 1, 1, 1, 1, 1,): [28, -14, 10, -6, 4, -2, 3, -1, 0, 2, 1, 1, -2, 0,
                                          0, 0, 0, 2, -1, 1, 0, -4, -1, 1, -1, 0, 0, -1, 2, 1],
                 (2, 2, 1, 1, 1, 1, 1,): [27, -15, 9, -5, 7, -3, 2, 0, 1, -3, 0, 0, 1, 1,
                                          -1, 1, -1, -1, 0, 0, -1, 3, 0, 0, 0, 1, -1, 2, -3, 0],
                 (2, 1, 1, 1, 1, 1, 1, 1,): [8, -6, 5, -4, 4, -3, 3, -2, 2, -2, 2, -1, 1, -1,
                                             1] + [0] * 7 + [1, -1] * 4,
                 (1, 1, 1, 1, 1, 1, 1, 1, 1,): [1, -1] * 10 + [1] + [1, -1] * 4 + [1]}

        self.character_1 = _ch_1
        self.character_2 = _ch_2
        self.character_3 = _ch_3
        self.character_4 = _ch_4
        self.character_5 = _ch_5
        self.character_6 = _ch_6
        self.character_7 = _ch_7
        self.character_8 = _ch_8
        self.character_9 = _ch_9
        self.max_s_n = 9  # 因为就录入到9

    def get_s_n_character(self, s_n):
        if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
            err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
            logger.error(err_msg)
            return False
        if s_n > self.max_s_n:
            err_msg = "s_n={} must <= self.max_s_n={}".format(s_n, self.max_s_n)
            logger.error(err_msg)
            return False
        rst = eval("self.character_{}".format(s_n))
        return copy.deepcopy(rst)


class GiData(object):
    """Littlewood的书中称Order，本论文以及本证函数法中order另有所指"""

    def __init__(self):
        self.gi_1 = [1]
        self.gi_2 = [1, 1]
        self.gi_3 = [1, 3, 2]
        self.gi_4 = [1, 6, 8, 6, 3]
        self.gi_5 = [1, 10, 20, 30, 15, 20, 24]
        self.gi_6 = [1, 15, 40, 90, 45, 120, 144, 120, 90, 15, 40]
        self.gi_7 = [1, 21, 70, 210, 105, 420, 504, 840, 630, 105, 280, 504, 210, 420, 720]
        self.gi_8 = [1, 28, 112, 420, 210, 1120, 1344, 3360, 2520, 420, 1120, 4032, 1680, 3360, 5760, 5040, 1260,
                     1260, 3360, 1120, 2688, 105]
        self.gi_9 = [1, 36, 168, 756, 378, 2520, 3024, 10080, 7560, 1260, 3360, 18144, 7560, 15120, 25920, 45360,
                     11340, 11340, 30240, 10080, 24192, 945, 18144, 2240, 20160, 15120, 25920, 9072, 2520, 40320]
        self.max_s_n = 9  # 因为就录入到9

    def get_s_n_gi(self, s_n):
        if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
            err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
            logger.error(err_msg)
            return False
        if s_n > self.max_s_n:
            err_msg = "s_n={} must <= self.max_s_n={}".format(s_n, self.max_s_n)
            logger.error(err_msg)
            return False
        rst = eval("self.gi_{}".format(s_n))
        return copy.deepcopy(rst)


def create_characters_and_gi(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储特征标矩阵和已经对齐的gi实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_characters_and_gi get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_characters_and_gi_finish_s_n()
    if not flag:
        err_msg = "get characters_and_gi finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} characters_and_gi had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc characters_and_gi s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的杨图并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        flag, yd_list_s_i = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db for create_characters_and_gi meet error with s_n={}, msg={}".format(
                s_n, yd_list_s_i)
            logger.error(err_msg)
            return False, err_msg
        flag, s_i_characters_and_gi = calc_single_characters_and_gi(s_i, yd_list=yd_list_s_i)
        s_i_speed_time = int(time.time() - s_i_start_time)
        if not flag:
            err_msg = "calc s_i_characters_and_gi meet error with s_i={}, msg={}".format(s_i, s_i_characters_and_gi)
            logger.error(err_msg)
            return False, err_msg
        if not s_i_characters_and_gi:
            err_msg = "calc s_i_characters_and_gi should not kill by recursion with s_i={}, pls check".format(s_i)
            logger.error(err_msg)
            return False, err_msg
        # 既然没问题了，那就存（别忘了也要更新Finish_Sn）
        flag, msg = save_characters_and_gi(s_i, s_i_characters_and_gi, yd_list_s_i, s_i_speed_time)
        if not flag:
            err_msg = "save s_i_characters_and_gi meet error with s_i={}, msg={}".format(s_i, s_i_characters_and_gi)
            logger.error(err_msg)
            return False, err_msg
        flag, msg = save_characters_and_gi_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_characters_and_gi_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_characters_and_gi s_n from {} to {} done, return True, finish_s_n={}, "
                "using time={}s".format(finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_characters_and_gi(s_n: int, characters_and_gi_dict: dict, yd_list: list, speed_time: int):
    """
    它们的落盘格式为：
    <CG>/characters_and_gi_info/Sn.pkl       ->
    {
    "file_name": "Sn",
    "data": {"character": character_matrix, "gi": gi_array},
    "flags": {"speed_time": speed_time,
              "young_diagram_index": young_diagram_list_by_yamanouchi}  # 作为参考，非数据
    }

    其中，
    Sn表示n阶置换群;
    character_matrix是按照Littlewood书中给的表格和Yamanouchi序存放特征标矩阵;  # np.ndarray(int)
    gi是已经和特征标矩阵的列对齐的gi，同样来自Littlewood的书中;                 # np.ndarray(int)
    speed_time表示计算用时（秒）

    例如：
    <CG>/characters_and_gi_info/S4.pkl:
    {"character": [[ 1  1  1  1  1]
                   [ 3  1  0 -1 -1]
                   [ 2  0 -1  0  2]
                   [ 3 -1  0  1 -1]
                   [ 1 -1  1 -1  1]],
     "gi": [1, 6, 8, 6, 3]}
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(characters_and_gi_dict, dict):
        err_msg = "characters_and_gi_dict={} with type={} must be list".format(
            characters_and_gi_dict, type(characters_and_gi_dict))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd_list, list):
        err_msg = "get yd_list={} with type={} must be list".format(yd_list, type(yd_list))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = CharacterAndGiInfo(s_n)
    _, file_name = get_characters_and_gi_file_name(s_n)
    table = {"file_name": file_name,
             "data": characters_and_gi_dict,
             "flags": {"speed_time": speed_time,
                       "young_diagram_index": yd_list}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_characters_and_gi_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_characters_and_gi_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = CharacterAndGiInfo(0)
    _, finish_file_name = get_characters_and_gi_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {
                           "S{}".format(s_n): s_n_speed_time
                       }}}

    if finish_s_n_before == min_s_n_of_characters_and_gi - 1:
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table, point_key="flags")
        if not flag:
            return flag, msg
    else:
        flag, data = db_info.query_by_file_name(finish_file_name)
        if not flag:
            return flag, data
        history_times = data.get("flags", {}).get("history_times")
        if not history_times or not isinstance(history_times, dict):  # 有，就不应该是空字典
            err_msg = "old history_times={} must real dict, but not, with data={}".format(history_times, data)
            logger.error(err_msg)
            return False, err_msg
        # 更新table里的历史时间
        history_times.update({"S{}".format(s_n): s_n_speed_time})
        table["flags"]["history_times"] = history_times
        flag, msg = db_info.update_by_file_name(finish_file_name, partial_table=table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(finish_file_name, partial_table=table, point_key="flags")
        if not flag:
            return flag, msg

    return True, None


def get_characters_and_gi_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 1
    """
    _, finish_file_name = get_characters_and_gi_finish_s_n_name()
    flag, data = CharacterAndGiInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get characters_and_gi finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_characters_and_gi - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_characters_and_gi:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_characters_and_gi,
                                                                          data)
        return False, err_msg


def load_characters_and_gi(s_n: int, is_flag_true_if_not_s_n=True):
    """
    取得s_n下的特征标矩阵和已经对齐的gi
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_characters_and_gi_file_name(s_n)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = CharacterAndGiInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query characters_and_gi with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        characters_and_gi_dict = data.get("data")
        if characters_and_gi_dict:
            # 只检查有没有 不对内容做检查了
            return True, characters_and_gi_dict  # bingo！

        else:
            err_msg = "characters_and_gi queried from db, but cannot get characters_and_gi from data" \
                      "with data={}, characters_and_gi={} from db".format(data, characters_and_gi_dict)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist characters_and_gi db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def calc_single_characters_and_gi(s_n: int, yd_list: list=None):
    """
    按照Yamanouchi序把录入的字典转换成矩阵，同时取得gi
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_characters_and_gi:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_characters_and_gi)
        logger.error(err_msg)
        return False, err_msg
    if yd_list is None:
        flag, yd_list = load_young_diagrams(s_n, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db for single_characters_and_gi meet error with s_n={}, msg={}".format(
                s_n, yd_list)
            logger.error(err_msg)
            return False, err_msg
    else:  # 如果是获得的yd_list，采取信任态度
        if not isinstance(yd_list, list):
            err_msg = "get yd_list={} with type={} must be list".format(yd_list, type(yd_list))
            logger.error(err_msg)
            return False, err_msg

    character_data = CharacterData()
    gi_data = GiData()
    if s_n > character_data.max_s_n or s_n > gi_data.max_s_n:
        err_msg = "s_n={} must <= character_data.max_s_n={} and <= gi_data.max_s_n={}".format(
            s_n, character_data.max_s_n, gi_data.max_s_n)
        logger.error(err_msg)
        return False, err_msg

    character_dict = character_data.get_s_n_character(s_n)
    gi_list = gi_data.get_s_n_gi(s_n)
    if not character_dict or not gi_list:
        err_msg = "cannot get character_dict={} and gi_list={}".format(character_dict, gi_list)
        logger.error(err_msg)
        return False, err_msg

    # 组装
    character_list = [character_dict[tuple(i)] for i in yd_list]
    character_matrix = np.array(character_list, dtype=int)
    gi_array = np.array(gi_list)
    characters_and_gi_dict = {"character": character_matrix,
                              "gi": gi_array}

    return True, characters_and_gi_dict
