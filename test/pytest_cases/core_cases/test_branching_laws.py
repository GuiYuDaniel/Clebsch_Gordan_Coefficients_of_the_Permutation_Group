# -*- coding:utf8 -*-
"""
测试
core/branching_laws.py
下所有功能是否正确执行
"""


import os
import pytest
import time
from conf.cgc_config import default_s_n, cgc_rst_folder
# from core.cgc_utils.cgc_local_db import get_branching_laws_file_name, get_branching_laws_finish_s_n_name
from core.branching_laws import calc_single_branching_law, create_branching_laws
# from core.branching_laws import load_branching_laws, get_branching_laws_finish_s_n
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestBranchingLaws(object):

    def setup_class(self):
        # regular BL
        self.yd_1 = [2, 1]  # single young_diagram
        self.bl_1 = (2, [1, 0], [0, 1], [[2], [1, 1]],)
        self.yd_2 = [8, 4, 4, 1]
        self.bl_2 = (3, [3, 2, 0], [0, 3, 7], [[8, 4, 4], [8, 4, 3, 1], [7, 4, 4, 1]])
        self.yd_3 = [1]
        self.bl_3 = (1, [0], [0], [[]],)  # 唯一一个特殊的before_YD
        # not regular BL
        self.yd_n1 = [2, 5]
        self.bl_n1 = (1, [1], [4], [[2, 4]])
        self.yd_n2 = [1, 2, 1]
        self.bl_n2 = (2, [2, 1], [0, 1], [[1, 2], [1, 1, 1]])
        # error BL
        self.yd_e1 = []
        self.yd_e2 = [[2, 1]]
        self.yd_e3 = {0: 2, 1: 1}
        self.yd_e4 = [4, "2", 1]
        self.yd_e5 = [-2, -3]
        self.yd_e6 = [2, 1.1]
        self.yd_e7 = [1, 0]

    def teardown_class(self):
        pass

    def test_calc_single_branching_law(self):
        # regular
        for i in [1, 2, 3]:
            yd = eval("self.yd_{}".format(i))
            bl = eval("self.bl_{}".format(i))
            flag, data = calc_single_branching_law(yd)
            assert flag
            assert data == bl

        # not regular
        for i in [1, 2]:
            yd = eval("self.yd_n{}".format(i))
            bl = eval("self.bl_n{}".format(i))
            flag, data = calc_single_branching_law(yd)
            assert flag
            assert data == bl

        # must error
        for i in range(1, 8):
            yd = eval("self.yd_e{}".format(i))
            flag, data = calc_single_branching_law(yd)
            logger.info("Error is supposed here!")
            assert not flag
            assert isinstance(data, str)
