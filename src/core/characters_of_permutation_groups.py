# -*- coding: utf-8 -*-
"""
this code for creating Characters of Permutation Groups by Sn
计算Sn的置换群特征标
"""

# 就是表示矩阵对角线上元素之和（trace迹）。对于有限群，双不表示的特征标正交归g
# TODO 特征标表很有意思的，意义也要再看看https://blog.csdn.net/wwxy1995/article/details/103744747


import time
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name, get_young_diagrams_finish_s_n_name
from conf.cgc_config import default_s_n
from utils.log import get_logger


logger = get_logger(__name__)


def func_pass():
    # TODO 先跳过，不一定真算，也可以拿来主义，反正结果有很多
    pass
