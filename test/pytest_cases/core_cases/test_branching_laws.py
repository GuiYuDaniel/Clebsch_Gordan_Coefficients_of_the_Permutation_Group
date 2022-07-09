# -*- coding:utf8 -*-
"""
测试
core/branching_laws.py
下所有功能是否正确执行
"""


import copy
import pytest
from conf.cgc_config import cgc_rst_folder
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo, BranchingLawInfo
from core.cgc_utils.cgc_local_db import get_branching_laws_file_name, get_branching_laws_finish_s_n_name
from core.branching_laws import calc_single_branching_law, create_branching_laws
from core.branching_laws import load_branching_law, get_branching_laws_finish_s_n
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestBranchingLaws(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_branching_laws_protected")
        self.protector.protector_setup()

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

        # 准备young_diagrams
        flag, msg = create_young_diagrams(4)
        assert flag
        assert msg == 4
        _, self.s_1_yd_list = load_young_diagrams(1, is_flag_true_if_not_s_n=True)
        _, self.s_2_yd_list = load_young_diagrams(2, is_flag_true_if_not_s_n=True)
        _, self.s_3_yd_list = load_young_diagrams(3, is_flag_true_if_not_s_n=True)
        _, self.s_4_yd_list = load_young_diagrams(4, is_flag_true_if_not_s_n=True)
        self.create_time_dict = {}  # 用于检查计算好的部分不会重复计算

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    def test_001_create_branching_laws_s_n_1_to_2(self):
        # check with no db
        flag, bl = load_branching_law(1, self.yd_3, is_return_tuple=True)
        assert flag
        assert bl is False
        flag, finish_s_n = get_branching_laws_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # create
        flag, data = create_branching_laws(2)
        assert flag
        assert data == 2

        # check answers
        flag, bl = load_branching_law(1, self.yd_3, is_return_tuple=True)
        assert flag
        assert bl == self.bl_3
        flag, bl = load_branching_law(3, self.yd_1, is_return_tuple=True)
        assert flag
        assert bl is False
        _, file_name = get_branching_laws_file_name(1, [1])
        flag, data = BranchingLawInfo(1).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        self.create_time_dict[file_name] = data.get("create_time")

        # check finish s_n
        flag, finish_s_n = get_branching_laws_finish_s_n()
        assert flag
        assert finish_s_n == 2

        # history_times
        _, finish_file_name = get_branching_laws_finish_s_n_name()
        flag, data = BranchingLawInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YoungDiagramInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_002_create_branching_laws_s_n_3_to_4(self):
        # create
        flag, data = create_branching_laws(4)
        assert flag
        assert data == 4

        # check answers
        flag, bl = load_branching_law(1, self.yd_3, is_return_tuple=True)
        assert flag
        assert bl == self.bl_3
        flag, bl = load_branching_law(3, self.yd_1, is_return_tuple=True)
        assert flag
        assert bl == self.bl_1
        _, file_name = get_branching_laws_file_name(1, [1])
        flag, data = BranchingLawInfo(1).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get(file_name) == data.get("create_time")

        # check finish s_n
        flag, finish_s_n = get_branching_laws_finish_s_n()
        assert flag
        assert finish_s_n == 4

        # history_times
        _, finish_file_name = get_branching_laws_finish_s_n_name()
        flag, data = BranchingLawInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YoungDiagramInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_calc_single_branching_law(self):
        # regular
        for i in [1, 2, 3]:
            yd = eval("self.yd_{}".format(i))
            yd_copy = copy.deepcopy(yd)
            bl = eval("self.bl_{}".format(i))
            flag, data = calc_single_branching_law(yd)
            assert flag
            assert data == bl
            assert yd == yd_copy  # 保证入参不变

        # not regular
        for i in [1, 2]:
            yd = eval("self.yd_n{}".format(i))
            yd_copy = copy.deepcopy(yd)
            bl = eval("self.bl_n{}".format(i))
            flag, data = calc_single_branching_law(yd)
            assert flag
            assert data == bl
            assert yd == yd_copy  # 保证入参不变

        # must error
        for i in range(1, 8):
            yd = eval("self.yd_e{}".format(i))
            yd_copy = copy.deepcopy(yd)
            flag, data = calc_single_branching_law(yd)
            logger.info("Error is supposed here!")
            assert not flag
            assert isinstance(data, str)
            assert yd == yd_copy  # 保证入参不变
