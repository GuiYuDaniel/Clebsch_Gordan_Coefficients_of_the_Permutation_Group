# -*- coding:utf8 -*-
"""
测试
core/young_tableaux.py
下所有功能是否正确执行
"""


import copy
import pytest
from conf.cgc_config import cgc_db_name
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo, BranchingLawInfo, YoungTableInfo
from core.cgc_utils.cgc_local_db import get_young_tableaux_file_name, get_young_tableaux_finish_s_n_name
from core.branching_laws import create_branching_laws
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from core.young_tableaux import create_young_tableaux, calc_single_young_table
from core.young_tableaux import load_young_table, get_young_tableaux_finish_s_n
from core.young_tableaux import load_young_table_num, load_young_table_phase_factor
from core.young_tableaux import read_young_table_in_decreasing_page_order
from core.young_tableaux import quickly_calc_young_table_in_decreasing_page_order
from core.young_tableaux import is_young_table
from core.young_tableaux import calc_young_diagram_from_young_table, calc_s_n_from_young_table
from core.young_tableaux import get_s_i_index_in_young_table
from core.young_tableaux import _calc_single_phase_factor, calc_single_yt_phase_factor
from core.young_tableaux import quick_calc_phase_factor_list
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


def normal_calc_phase_factor_list(s_n: int, young_tableaux: dict):
    """平庸的逐个循环算相位因子Λ"""
    # 保护
    if s_n == 1:
        return True, [1]

    # 开始计算
    phase_factor_list = []
    for i in range(len(young_tableaux)):
        m = str(i + 1)
        yt = young_tableaux[m]
        yt_phase_factor = calc_single_yt_phase_factor(yt, _is_check_yt=False)
        phase_factor_list.append(yt_phase_factor)
    return True, phase_factor_list


class TestYoungTableaux(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_db_name, extension_name=".test_young_tableaux_protected")
        self.protector.protector_setup()

        # regular YT
        self.sn_1, self.yd_1 = 3, [2, 1]
        self.yt_1 = {"1": [[1, 2], [3]], "2": [[1, 3], [2]]}
        self.yt_1_Λ = [1, -1]
        self.sn_2, self.yd_2 = 1, [1]
        self.yt_2 = {"1": [[1]]}
        self.yt_2_Λ = [1]
        self.sn_3, self.yd_3 = 4, [3, 1]
        self.yt_3 = {"1": [[1, 2, 3], [4]],
                     "2": [[1, 2, 4], [3]],
                     "3": [[1, 3, 4], [2]]}
        self.yt_3_Λ = [1, -1, 1]
        # error YT
        self.sn_e1, self.yd_e1 = 3, [2, 2]

        # other YT  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-1
        self.regular_young_table_and_idp_order_dict = {
            "[[1, 2], [3]]": 1,
            "[[1, 3], [2]]": 2,
            "[[1, 2, 3], [4]]": 1,
            "[[1, 2, 4], [3]]": 2,
            "[[1, 3, 4], [2]]": 3,
            "[[1, 2], [3, 4]]": 1,
            "[[1, 3], [2, 4]]": 2,
            "[[1, 3], [2], [4]]": 2,
            "[[1, 2, 4, 5], [3]]": 3,
            "[[1, 2, 4], [3, 5]]": 2,
            "[[1, 2, 5], [3, 4]]": 4,
            "[[1, 4, 5], [2], [3]]": 6,
            "[[1, 4], [2, 5], [3]]": 5,
            "[[1, 4], [2], [3], [5]]": 3,
            "[[1, 2, 4, 5, 6], [3]]": 4,
            "[[1, 3, 5, 6], [2, 4]]": 9,
            "[[1, 2, 3, 5], [4], [6]]": 2,
            "[[1, 3, 4, 6], [2], [5]]": 7,
            "[[1, 2, 3], [4, 5, 6]]": 1,
            "[[1, 3, 5], [2, 4, 6]]": 5,
            "[[1, 3, 5], [2, 4], [6]]": 5,
            "[[1, 2, 4], [3, 6], [5]]": 7,
            "[[1, 4, 5], [2, 6], [3]]": 11,
            "[[1, 4], [2, 5], [3, 6]]": 5,
            "[[1, 2, 3], [4], [5], [6]]": 1,
            "[[1, 3], [2, 5], [4], [6]]": 4,
            "[[1, 3], [2, 6], [4], [5]]": 7,
            "[[1, 6], [2], [3], [4], [5]]": 5,
        }

        # other_wrong_yt
        self.wrong_yts = (None, ([1, 2], [3],), [], [1, 2], [(1, 3,), [2]], [None, None], [["1"], ["2"]],
                          [[2, 1]], [[1, 3], [4]], [[-3, -2], [-1]])

        # 准备前文
        flag, msg = create_young_diagrams(4)
        assert flag
        assert msg == 4
        flag, msg = create_branching_laws(4)
        assert flag
        assert msg == 4

        self.create_time_dict = {}

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    def test_001_create_young_tableaux_s_n_1_to_2(self):
        # check with no db
        for i in [1, 2, 3]:
            s_i, yd_i = eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i))
            flag, yt = load_young_table(s_i, yd_i)
            assert flag
            assert yt is False
        flag, finish_s_n = get_young_tableaux_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # create
        flag, data = create_young_tableaux(2)
        assert flag
        assert data == 2

        # check answers
        for i in [1, 2, 3]:
            s_i, yd_i, yt_i = eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i)), eval("self.yt_{}".format(i))
            yt_i_Λ = eval("self.yt_{}_Λ".format(i))
            flag, yt = load_young_table(s_i, yd_i)
            assert flag
            if s_i > 2:
                assert yt is False
            else:
                assert yt == yt_i
                _, file_name = get_young_tableaux_file_name(s_i, yd_i)
                flag, data = YoungTableInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == yt_i
                assert data.get("flags").get("speed_time") >= 0
                flag, yt_num = load_young_table_num(s_i, yd_i)
                assert flag
                assert isinstance(yt_num, int)
                assert yt_num == len(yt_i) == data.get("flags").get("total_num")
                flag, yt_Λ = load_young_table_phase_factor(s_i, yd_i)
                assert flag
                assert isinstance(yt_Λ, list)
                assert yt_Λ == yt_i_Λ
                self.create_time_dict[file_name] = data.get("create_time")

        # check finish s_n
        flag, finish_s_n = get_young_tableaux_finish_s_n()
        assert flag
        assert finish_s_n == 2

        # history_times
        _, finish_file_name = get_young_tableaux_finish_s_n_name()
        flag, data = YoungTableInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YoungTableInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_002_create_young_tableaux_s_n_3_to_4(self):
        # create
        flag, data = create_young_tableaux(4)
        assert flag
        assert data == 4

        # check answers
        for i in [1, 2, 3]:
            s_i, yd_i, yt_i = eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i)), eval("self.yt_{}".format(i))
            yt_i_Λ = eval("self.yt_{}_Λ".format(i))
            flag, yt = load_young_table(s_i, yd_i)
            assert flag
            if s_i > 4:
                assert yt is False
            else:
                assert yt == yt_i
                _, file_name = get_young_tableaux_file_name(s_i, yd_i)
                flag, data = YoungTableInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == yt_i
                assert data.get("flags").get("speed_time") >= 0
                flag, yt_num = load_young_table_num(s_i, yd_i)
                assert flag
                assert isinstance(yt_num, int)
                assert yt_num == len(yt_i) == data.get("flags").get("total_num")
                flag, yt_Λ = load_young_table_phase_factor(s_i, yd_i)
                assert flag
                assert isinstance(yt_Λ, list)
                assert yt_Λ == yt_i_Λ
                if s_i <= 2:
                    assert self.create_time_dict.get(file_name) == data.get("create_time")
                else:
                    self.create_time_dict[file_name] = data.get("create_time")

        # check finish s_n
        flag, finish_s_n = get_young_tableaux_finish_s_n()
        assert flag
        assert finish_s_n == 4

        # history_times
        _, finish_file_name = get_young_tableaux_finish_s_n_name()
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

    def test_calc_single_young_table(self):
        # regular
        for i in [1, 2, 3]:
            s_i, yd_i, yt_i = eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i)), eval("self.yt_{}".format(i))
            yt_i_copy = copy.deepcopy(yt_i)
            flag, data = calc_single_young_table(s_i, yd_i)
            assert flag
            assert data == yt_i
            assert yt_i == yt_i_copy  # 保证入参不变

        # regular
        for i in [1]:
            s_ei, yd_ei = eval("self.sn_e{}".format(i)), eval("self.yd_e{}".format(i))
            flag, data = calc_single_young_table(s_ei, yd_ei)
            logger.info("Error is supposed here!")
            assert not flag
            assert isinstance(data, str)

    def test_is_young_table(self):
        # 正例1
        head_s_n = 1
        tail_s_n = 3
        for s_n in range(head_s_n, tail_s_n + 1):
            yt_dict = eval("self.yt_{}".format(s_n))
            for idp_order_str, yt in yt_dict.items():
                assert is_young_table(yt)

        # 正例2
        for yt_str, idp_order in self.regular_young_table_and_idp_order_dict.items():
            assert is_young_table(eval(yt_str))

        # 反例
        for wrong_yt in self.wrong_yts:
            assert not is_young_table(wrong_yt)
            logger.info("Error is supposed here!")

    def test_calc_young_diagram_from_young_table(self):
        # 正
        head_s_n = 1
        tail_s_n = 3
        for s_n in range(head_s_n, tail_s_n + 1):
            yt_dict = eval("self.yt_{}".format(s_n))
            yd_answer = eval("self.yd_{}".format(s_n))
            for idp_order_str, yt in yt_dict.items():
                flag, yd = calc_young_diagram_from_young_table(yt)
                assert flag
                assert yd_answer == yd

        # 反
        for wrong_yt in self.wrong_yts:
            flag, yd = calc_young_diagram_from_young_table(wrong_yt)
            logger.info("Error is supposed here!")
            assert not flag
            assert yd is None

    def test_calc_s_n_from_young_table(self):
        # 正
        head_s_n = 1
        tail_s_n = 3
        for s_n in range(head_s_n, tail_s_n + 1):
            yt_dict = eval("self.yt_{}".format(s_n))
            s_n_answer = eval("self.sn_{}".format(s_n))
            for idp_order_str, yt in yt_dict.items():
                flag, s_n_calc = calc_s_n_from_young_table(yt)
                assert flag
                assert s_n_answer == s_n_calc

        # 反
        for wrong_yt in self.wrong_yts:
            flag, s_n_calc = calc_s_n_from_young_table(wrong_yt)
            logger.info("Error is supposed here!")
            assert not flag
            assert s_n_calc is None

    def test_get_s_i_index_in_young_table(self):
        # 原函数
        # 寻找矩阵c中元素d坐标的函数
        def FindSeat(Matrix, value):
            """
            寻找矩阵c中元素d坐标的函数
            :param c:Mateix
            :param d: value
            :return: Seat of value
            """
            x = 0
            y = 0
            for i in Matrix:
                try:
                    y = i.index(value)
                except:
                    x += 1
                else:
                    break
            Seat = [x, y]
            return Seat
        # 格式：yt, point, row_answer, col_answer
        data_1 = ([[1, 2], [3]], 1, 0, 0)
        data_2 = ([[1, 2], [3]], 2, 0, 1)
        data_3 = ([[1, 2], [3]], 3, 1, 0)
        data_4 = ([[1, 2], [3, 4]], 4, 1, 1)
        for i in range(1, 5):  # [1, 2, 3, 4]
            yt, point, row_answer, col_answer = eval("data_{}".format(i))
            row_1, col_1 = FindSeat(yt, point)
            row_2, col_2 = get_s_i_index_in_young_table(point, yt)
            assert row_1 == row_2 == row_answer
            assert col_1 == col_2 == col_answer

    def test_read_young_table_in_decreasing_page_order(self):
        # 正
        for idp_order_str, yt in self.yt_3.items():
            flag, rst = read_young_table_in_decreasing_page_order(yt, self.yt_3)
            assert flag
            assert rst == idp_order_str

        # 反
        flag, rst = read_young_table_in_decreasing_page_order(self.yt_2.get("1"), self.yt_3)
        assert flag
        assert rst is None

    def test_003_quickly_calc_young_table_in_decreasing_page_order(self):
        # 准备前文
        flag, msg = create_young_diagrams(6)
        assert flag
        assert msg == 6
        flag, msg = create_branching_laws(6)
        assert flag
        assert msg == 6
        flag, data = create_young_tableaux(6)
        assert flag
        assert msg == 6
        # flag, msg = create_young_diagrams(9)
        # assert flag
        # assert msg == 9
        # flag, msg = create_branching_laws(9)
        # assert flag
        # assert msg == 9
        # flag, data = create_young_tableaux(9)
        # assert flag
        # assert msg == 9

        for yt_str, idp_order in self.regular_young_table_and_idp_order_dict.items():
            yt = eval(yt_str)
            flag, data = quickly_calc_young_table_in_decreasing_page_order(yt, is_rst_int=True)
            assert flag
            assert data == idp_order

        for idp_order_str, yt in self.yt_3.items():
            flag, rst = quickly_calc_young_table_in_decreasing_page_order(yt, yd=[3, 1], s_n=4)
            assert flag
            assert rst == idp_order_str

    def test__calc_single_phase_factor(self):
        # answer
        y_target_1 = [1, 3, 4, 6, 2, 5]
        phase_factor_1 = 1
        exchange_list_1 = [(4, 5), (3, 4), (2, 3), (5, 6)]

        y_target_2 = [1]
        phase_factor_2 = 1
        exchange_list_2 = []

        y_target_3 = [1, 4, 2, 5, 3, 6]
        phase_factor_3 = -1
        exchange_list_3 = [(2, 3), (4, 5), (3, 4)]
        
        for num in [1, 2, 3]:
            y_target = eval("y_target_{}".format(num))
            phase_factor = eval("phase_factor_{}".format(num))
            exchange_list = eval("exchange_list_{}".format(num))

            y_start = list(range(1, len(y_target)+1))
            assert _calc_single_phase_factor(y_target, y_start) == phase_factor
            assert y_target == y_start  # y_start会变

            y_start = list(range(1, len(y_target) + 1))
            phase_factor_rst, exchange_list_rst = _calc_single_phase_factor(y_target, y_start,
                                                                            is_output_exchange_list=True)
            assert phase_factor_rst == phase_factor
            assert exchange_list == exchange_list_rst
            assert y_target == y_start  # y_start会变

    def test_calc_single_yt_phase_factor(self):
        # answer
        yt_target_1 = [[1, 3, 4, 6],
                       [2, 5]]
        phase_factor_1 = 1
        exchange_list_1 = [(4, 5), (3, 4), (2, 3), (5, 6)]

        yt_target_2 = [[1]]
        phase_factor_2 = 1
        exchange_list_2 = []

        yt_target_3 = [[1, 4],
                       [2, 5],
                       [3],
                       [6]]
        phase_factor_3 = -1
        exchange_list_3 = [(2, 3), (4, 5), (3, 4)]

        for num in [1, 2, 3]:
            yt_target = eval("yt_target_{}".format(num))
            phase_factor = eval("phase_factor_{}".format(num))
            exchange_list = eval("exchange_list_{}".format(num))

            assert calc_single_yt_phase_factor(yt_target) == phase_factor
            assert calc_single_yt_phase_factor(yt_target, _is_check_yt=False) == phase_factor
            phase_factor_rst, exchange_list_rst = calc_single_yt_phase_factor(yt_target,
                                                                              is_output_exchange_list=True)
            assert phase_factor_rst == phase_factor
            assert exchange_list_rst == exchange_list
            phase_factor_rst, exchange_list_rst = calc_single_yt_phase_factor(yt_target, _is_check_yt=False,
                                                                              is_output_exchange_list=True)
            assert phase_factor_rst == phase_factor
            assert exchange_list_rst == exchange_list

    def test_004_quick_calc_phase_factor_list(self):
        for s_i in range(1, 6+1):
            _, young_diagrams = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
            for yd in young_diagrams:
                # 拿quick_calc_phase_factor_list的结果
                flag, phase_factor_list_quick = load_young_table_phase_factor(s_i, yd, is_flag_true_if_not_s_n=False)
                assert flag, "msg={}".format(phase_factor_list_quick)
                # 按照normal_calc_phase_factor_list算
                _, yt = load_young_table(s_i, yd, is_flag_true_if_not_s_n=False)
                assert flag
                flag, phase_factor_list_normal = normal_calc_phase_factor_list(s_i, yt)
                assert flag
                assert phase_factor_list_quick == phase_factor_list_normal

