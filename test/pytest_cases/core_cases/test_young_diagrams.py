# -*- coding:utf8 -*-
"""
测试
core/young_diagrams.py
下所有功能是否正确执行
"""


import os
import pytest
import time
from conf.cgc_config import default_s_n, cgc_rst_folder
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name, get_young_diagrams_finish_s_n_name
from core.young_diagrams import calc_single_young_diagrams, create_young_diagrams
from core.young_diagrams import load_young_diagrams, get_young_diagrams_finish_s_n
from core.young_diagrams import is_young_diagram, calc_s_n_from_young_diagram
from core.young_diagrams import calc_young_diagram_dagger
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestYoungDiagrams(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_young_diagrams_protected")
        self.protector.protector_setup()

        self.young_diagrams_s_1 = [[1]]
        self.young_diagrams_s_2 = [[2], [1, 1]]
        self.young_diagrams_s_3 = [[3], [2, 1], [1, 1, 1]]
        self.young_diagrams_s_4 = [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        self.young_diagrams_s_5 = [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        self.young_diagrams_s_6 = [[6], [5, 1],
                                   [4, 2], [4, 1, 1],
                                   [3, 3], [3, 2, 1], [3, 1, 1, 1],
                                   [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1, 1]]
        self.young_diagrams_s_7 = [[7], [6, 1],
                                   [5, 2], [5, 1, 1],
                                   [4, 3], [4, 2, 1], [4, 1, 1, 1],
                                   [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [3, 1, 1, 1, 1],
                                   [2, 2, 2, 1], [2, 2, 1, 1, 1], [2, 1, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1, 1, 1]]
        self.wrong_young_diagram_list = [None, [], [1, 2], [3, 0], [3, 2, -1], [-1, -2]]
        self.young_diagrams_dict = {"+": [[3, 2], [3, 1, 1], [2, 2], [1], [1, 1, 1, 1]],
                                    "-": [[2, 2, 1], [3, 1, 1], [2, 2], [1], [4]]}
        self.default_s_n = default_s_n
        _, self.s_1_file_name = get_young_diagrams_file_name(1)
        _, self.s_1_full_file_name = get_young_diagrams_file_name(1, is_full_path=True)
        _, self.s_n_finish_file_name = get_young_diagrams_finish_s_n_name()
        _, self.s_n_finish_full_file_name = get_young_diagrams_finish_s_n_name(is_full_path=True)

        _, self.s_4_full_file_name = get_young_diagrams_file_name(4, is_full_path=True)
        self.create_time_dict = {}  # 用于检查计算好的部分不会重复计算

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    def test_001_create_young_diagrams_s_n_1(self):  # calc and save, return True, None
        """for s_n=1, there is no finish db"""
        # check with no db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.s_1_full_file_name + ex)
            assert not os.path.exists(self.s_n_finish_full_file_name + ex)
        flag, young_diagrams = load_young_diagrams(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert young_diagrams is False
        flag, young_diagrams = load_young_diagrams(1, is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(young_diagrams, str)
        flag, finish_s_n = get_young_diagrams_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # check create_young_diagrams_s_n_1
        flag, msg = create_young_diagrams(1)
        assert flag
        assert msg == 1

        # check answer
        flag, data = YoungDiagramInfo(1).query_by_file_name(self.s_1_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        self.create_time_dict["S1"] = data.get("create_time")
        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_1_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_full_file_name + ex)
        flag, young_diagrams = load_young_diagrams(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert young_diagrams == self.young_diagrams_s_1

        # check finish s_n
        flag, finish_s_n = get_young_diagrams_finish_s_n()
        assert flag
        assert finish_s_n == 1

        # history_times
        _, finish_file_name = get_young_diagrams_finish_s_n_name()
        flag, data = YoungDiagramInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        assert isinstance(data.get("flags").get("history_times").get("S1"), int)
        assert 0 <= data.get("flags").get("history_times").get("S1") <= 1
        flag, data_txt = YoungDiagramInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        assert isinstance(data.get("history_times").get("S{}".format(1)), int)

    def test_002_create_young_diagrams_s_n_2_to_4(self):
        # check create_young_diagrams_s_n 2 to 4
        flag, msg = create_young_diagrams(4)
        assert flag
        assert msg == 4

        flag, data = YoungDiagramInfo(1).query_by_file_name(self.s_1_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get("S1") == data.get("create_time")

        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_4_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_full_file_name + ex)
        for i in [1, 2, 3, 4]:
            flag, young_diagrams = load_young_diagrams(i, is_flag_true_if_not_s_n=True)
            assert flag
            assert young_diagrams == eval("self.young_diagrams_s_{}".format(i))
        flag, finish_s_n = get_young_diagrams_finish_s_n()
        assert flag
        assert finish_s_n == 4
        _, finish_file_name = get_young_diagrams_finish_s_n_name()
        flag, data = YoungDiagramInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("flags").get("history_times").get("S{}".format(i)) is None
        flag, data_txt = YoungDiagramInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("history_times").get("S{}".format(i)) is None

    def test_003_calc_single_young_diagrams(self):
        """1 to 4 is finish, 5 to 7 is needing calc"""
        head_s_n = 1
        tail_s_n = 7
        head_time = time.time()
        for s_n in range(head_s_n, tail_s_n + 1):
            flag, young_diagrams = calc_single_young_diagrams(s_n, recursion_deep=7)
            logger.debug("for s_n={}, test young_diagrams={}".format(s_n, young_diagrams))
            assert flag
            assert young_diagrams == eval("self.young_diagrams_s_{}".format(s_n)), "s_n={}".format(s_n)
        tail_time = time.time()
        logger.debug("calc young_diagrams from S{} to S{}, used {}s".format(head_s_n, tail_s_n, tail_time - head_time))
        _, finish_file_name = get_young_diagrams_finish_s_n_name()
        flag, data = YoungDiagramInfo(0).query_by_file_name(finish_file_name)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("flags").get("history_times").get("S{}".format(i)) is None
        flag, data_txt = YoungDiagramInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("history_times").get("S{}".format(i)) is None

    def test_calc_single_young_diagrams_out_of_recursion_deep(self):
        s_n = 1000
        flag, young_diagrams = calc_single_young_diagrams(s_n, recursion_deep=7)
        logger.debug("for s_n={}, test young_diagrams={}".format(s_n, young_diagrams))
        assert flag
        assert young_diagrams is False

    def test_is_young_diagram(self):
        # 正例
        head_s_n = 1
        tail_s_n = 7
        for s_n in range(head_s_n, tail_s_n + 1):
            young_diagrams = eval("self.young_diagrams_s_{}".format(s_n))
            for yd in young_diagrams:
                assert is_young_diagram(yd)

        # 反例
        for wrong_yd in self.wrong_young_diagram_list:
            assert not is_young_diagram(wrong_yd)

    def test_calc_s_n_from_young_diagram(self):
        head_s_n = 1
        tail_s_n = 7
        for s_n in range(head_s_n, tail_s_n + 1):
            young_diagrams = eval("self.young_diagrams_s_{}".format(s_n))
            for yd in young_diagrams:
                flag, s_n_calc = calc_s_n_from_young_diagram(yd)
                assert flag
                assert s_n_calc == s_n

    def test_calc_young_diagram_dagger(self):
        yd_list_1, yd_list_2 = self.young_diagrams_dict["+"], self.young_diagrams_dict["-"]
        for yd_1, yd_2 in zip(yd_list_1, yd_list_2):
            flag, yd_1_dagger = calc_young_diagram_dagger(yd_1)
            assert flag
            assert yd_1_dagger == yd_2
            flag, yd_1_dagger_dagger = calc_young_diagram_dagger(yd_1_dagger)
            assert flag
            assert yd_1 == yd_1_dagger_dagger

        head_s_n = 1
        tail_s_n = 7
        for s_n in range(head_s_n, tail_s_n + 1):
            young_diagrams = eval("self.young_diagrams_s_{}".format(s_n))
            for yd in young_diagrams:
                flag, yd_dagger = calc_young_diagram_dagger(yd)
                assert flag
                assert is_young_diagram(yd_dagger)
                flag, yd_dagger_dagger = calc_young_diagram_dagger(yd_dagger)
                assert flag
                assert yd == yd_dagger_dagger

