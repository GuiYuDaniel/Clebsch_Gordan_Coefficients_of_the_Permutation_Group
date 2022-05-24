# -*- coding:utf8 -*-
"""
测试
core/eigenvalues.py
下所有功能是否正确执行
"""


import os
import time
import math
import numpy as np
from conf.cgc_config import cgc_rst_folder
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from core.characters_and_gi import create_characters_and_gi
from core.cg_order import create_cg_order, save_cg_order, load_cg_order
from core.cg_order import save_cg_order_finish_s_n, get_cg_order_finish_s_n
from core.cgc_utils.cgc_db_typing import CGOrderInfo
from core.cgc_utils.cgc_local_db import get_cg_order_file_name, get_cg_order_finish_s_n_name
from core.eigenvalues import calc_eigenvalues_of_2_cycle_class_of_s_n, _calc_single_eigenvalue_of_2_cycle_class_of_yd
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


# old function
# 二循环类算符对应杨盘的本征值
def geteigenvalue(r):
    result = 0
    lenth = len(r)
    for i in range(lenth):
        a = r[i] * (r[i]+1-2*(i+1)) / 2
        a = int(a)
        result += a
    return result


class TestEigenvalues(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_eigenvalues_protected")
        self.protector.protector_setup()

        # 准备前文
        flag, msg = create_young_diagrams(6)
        assert flag
        assert msg == 6

        # 数据来源见《群表示论的新途径》陈金全（上海科学技术出版社1984）图4.5-3
        # 格式 Sn, ev_list(对应young diagrams的 Yamanouchi序)
        self.ev_1 = [1]  # 注意，二循环类算符，本应从S2开始，其S1的值不是计算的，而是定义的延拓
        self.ev_2 = [1, -1]
        self.ev_3 = [3, 0, -3]
        self.ev_4 = [6, 2, 0, -2, -6]
        self.ev_5 = [10, 5, 2, 0, -2, -5, -10]
        self.ev_6 = [15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15]

        # 格式 Sn, ev, yd_list
        self.ev_to_yd_1 = (6, 3, [[4, 1, 1], [3, 3]])  # 注意，书中数据为了美观，这里没有按照Yamanouchi序绘图，我们按照了
        self.ev_to_yd_2 = (6, -3, [[3, 1, 1, 1], [2, 2, 2]])

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    def test__calc_single_eigenvalue_of_2_cycle_class_of_yd(self):
        for s_i in range(1, 6 + 1):
            ev_answer_list = eval("self.ev_{}".format(s_i))
            _, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
            assert len(ev_answer_list) == len(yd_list)
            for yd, ev_answer in zip(yd_list, ev_answer_list):
                flag, ev = _calc_single_eigenvalue_of_2_cycle_class_of_yd(yd)
                assert flag
                assert ev == ev_answer, "with s_i={}, yd={}".format(s_i, yd)
                if s_i > 1:
                    ev_old = geteigenvalue(yd)
                    assert ev == ev_old

    def test_calc_eigenvalues_of_2_cycle_class_of_s_n(self):
        for s_i in range(1, 6 + 1):
            ev_answer_list = eval("self.ev_{}".format(s_i))
            _, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
            assert len(ev_answer_list) == len(yd_list)
            flag, ev_list = calc_eigenvalues_of_2_cycle_class_of_s_n(yd_list)
            assert flag
            assert ev_list == ev_answer_list, "with s_i={}, yd_list={}".format(s_i, yd_list)

    def test_get_yd_list_by_eigenvalue(self):
        pass
