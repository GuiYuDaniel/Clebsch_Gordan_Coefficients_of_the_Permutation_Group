# -*- coding:utf8 -*-
"""
测试
core/yamanouchi_matrix.py
下所有功能是否正确执行
"""


import copy
import sympy as sp
import pytest
from conf.cgc_config import cgc_db_name
from core.cgc_utils.cgc_db_typing import YamanouchiMatrixInfo
from core.branching_laws import create_branching_laws
from core.young_diagrams import create_young_diagrams
from core.young_tableaux import create_young_tableaux
from core.yamanouchi_matrix import create_yamanouchi_matrix
from core.yamanouchi_matrix import get_yamanouchi_matrix_finish_s_n
from core.cgc_utils.cgc_local_db import get_yamanouchi_matrix_file_name, get_yamanouchi_matrix_finish_s_n_name
from core.yamanouchi_matrix import _calc_s_b_to_s_n_part
from core.branching_laws import load_branching_law
from core.young_tableaux import load_young_table, load_young_table_num
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestYamanouchiMatrixSP(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_db_name, extension_name=".test_yamanouchi_matrix_protected")
        self.protector.protector_setup()

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

        # (Sb, Sn)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
        self.ym_1 = sp.Matrix([[sp.Rational(-1)/2, sp.sqrt(3)/2], [sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_1, self.yd_1, self.ix_1 = 3, [2, 1], (2, 3,)
        self.ym_2 = sp.Matrix([[sp.Rational(-1)/3, sp.sqrt(8)/3, 0], [sp.sqrt(8)/3, sp.Rational(1)/3, 0], [0, 0, 1]])
        self.sn_2, self.yd_2, self.ix_2 = 4, [3, 1], (3, 4,)
        self.ym_3 = sp.Matrix([[-1, 0, 0], [0, sp.Rational(-1)/3, sp.sqrt(8)/3], [0, sp.sqrt(8)/3, sp.Rational(1)/3]])
        self.sn_3, self.yd_3, self.ix_3 = 4, [2, 1, 1], (3, 4,)
        self.ym_4 = sp.Matrix([[1, 0], [0, -1]])
        self.sn_4, self.yd_4, self.ix_4 = 4, [2, 2], (3, 4,)
        self.ym_5 = sp.Matrix([[sp.Rational(-1)/4, sp.sqrt(15)/4, 0, 0],
                               [sp.sqrt(15)/4, sp.Rational(1)/4, 0, 0],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])
        self.sn_5, self.yd_5, self.ix_5 = 5, [4, 1], (4, 5,)
        self.ym_6 = sp.Matrix([[1, 0, 0, 0, 0],
                               [0, sp.Rational(-1)/2, 0, sp.sqrt(3)/2, 0],
                               [0, 0, sp.Rational(-1)/2, 0, sp.sqrt(3)/2],
                               [0, sp.sqrt(3)/2, 0, sp.Rational(1)/2, 0],
                               [0, 0, sp.sqrt(3)/2, 0, sp.Rational(1)/2]])
        self.sn_6, self.yd_6, self.ix_6 = 5, [3, 2], (4, 5,)
        self.ym_7 = sp.Matrix([[-1, 0, 0, 0, 0, 0],
                               [0, sp.Rational(-1)/4, 0, sp.sqrt(15)/4, 0, 0],
                               [0, 0, -1/4, 0, sp.sqrt(15)/4, 0],
                               [0, sp.sqrt(15)/4, 0, sp.Rational(1)/4, 0, 0],
                               [0, 0, sp.sqrt(15)/4, 0, sp.Rational(1)/4, 0],
                               [0, 0, 0, 0, 0, 1]])
        self.sn_7, self.yd_7, self.ix_7 = 5, [3, 1, 1], (4, 5,)

        # (1, 2)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
        self.ym_11 = sp.Matrix([[1, 0], [0, -1]])
        self.sn_11, self.yd_11, self.ix_11 = 3, [2, 1], (1, 2,)
        self.ym_12 = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.sn_12, self.yd_12, self.ix_12 = 4, [3, 1], (1, 2,)
        self.ym_13 = sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
        self.sn_13, self.yd_13, self.ix_13 = 4, [2, 1, 1], (1, 2,)
        self.ym_14 = sp.Matrix([[1, 0], [0, -1]])
        self.sn_14, self.yd_14, self.ix_14 = 4, [2, 2], (1, 2,)
        self.ym_15 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
        self.sn_15, self.yd_15, self.ix_15 = 5, [4, 1], (1, 2,)
        self.ym_16 = sp.Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, -1]])
        self.sn_16, self.yd_16, self.ix_16 = 5, [3, 2], (1, 2,)
        self.ym_17 = sp.Matrix([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0],
                                [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]])
        self.sn_17, self.yd_17, self.ix_17 = 5, [3, 1, 1], (1, 2,)

        # (2, 3)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
        self.ym_21 = sp.Matrix([[sp.Rational(-1)/2, sp.sqrt(3)/2], [sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_21, self.yd_21, self.ix_21 = 3, [2, 1], (2, 3,)
        self.ym_22 = sp.Matrix([[1, 0, 0], [0, sp.Rational(-1)/2, sp.sqrt(3)/2], [0, sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_22, self.yd_22, self.ix_22 = 4, [3, 1], (2, 3,)
        self.ym_23 = sp.Matrix([[sp.Rational(-1)/2, sp.sqrt(3)/2, 0], [sp.sqrt(3)/2, sp.Rational(1)/2, 0], [0, 0, -1]])
        self.sn_23, self.yd_23, self.ix_23 = 4, [2, 1, 1], (2, 3,)
        self.ym_24 = sp.Matrix([[sp.Rational(-1)/2, sp.sqrt(3)/2], [sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_24, self.yd_24, self.ix_24 = 4, [2, 2], (2, 3,)
        self.ym_25 = sp.Matrix([[1, 0, 0, 0],
                                [0, 1, 0, 0],
                                [0, 0, sp.Rational(-1)/2, sp.sqrt(3)/2],
                                [0, 0, sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_25, self.yd_25, self.ix_25 = 5, [4, 1], (2, 3,)
        self.ym_26 = sp.Matrix([[1, 0, 0, 0, 0],
                                [0, sp.Rational(-1)/2, sp.sqrt(3)/2, 0, 0],
                                [0, sp.sqrt(3)/2, sp.Rational(1)/2, 0, 0],
                                [0, 0, 0, sp.Rational(-1)/2, sp.sqrt(3)/2],
                                [0, 0, 0, sp.sqrt(3)/2, sp.Rational(1)/2]])
        self.sn_26, self.yd_26, self.ix_26 = 5, [3, 2], (2, 3,)
        self.ym_27 = sp.Matrix([[1, 0, 0, 0, 0, 0],
                                [0, sp.Rational(-1)/2, sp.sqrt(3)/2, 0, 0, 0],
                                [0, sp.sqrt(3)/2, sp.Rational(1)/2, 0, 0, 0],
                                [0, 0, 0, sp.Rational(-1)/2, sp.sqrt(3)/2, 0],
                                [0, 0, 0, sp.sqrt(3)/2, sp.Rational(1)/2, 0],
                                [0, 0, 0, 0, 0, -1]])
        self.sn_27, self.yd_27, self.ix_27 = 5, [3, 1, 1], (2, 3,)
        pass

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    def test_001_create_yamanouchi_matrix_sn_2(self):
        # check with no db
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == 1

        # create
        flag, test_s_n = create_yamanouchi_matrix(2)
        assert flag
        assert test_s_n == 2

        # check answers
        all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
        for i in all_ij_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0
        all_in_num_list = list(range(1, 7 + 1))
        for i in all_in_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == 2

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [2]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [2]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_002_create_yamanouchi_matrix_sn_3_to_4(self):
        # create
        flag, test_s_n = create_yamanouchi_matrix(4)
        assert flag
        assert test_s_n == 4

        # check answers
        all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
        for i in all_ij_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0
        all_in_num_list = list(range(1, 7 + 1))
        for i in all_in_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == 4

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(2, test_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(2, test_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_002_create_yamanouchi_matrix_sn_5_to_6(self):
        # create
        flag, test_s_n = create_yamanouchi_matrix(6)
        assert flag
        assert test_s_n == 6

        # check answers
        all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
        for i in all_ij_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0
        all_in_num_list = list(range(1, 7 + 1))
        for i in all_in_num_list:
            ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
            ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
            if s_i <= test_s_n:
                _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                assert data.get("data") == ym
                assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == 6

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(2, test_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(2, test_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test__calc_s_b_to_s_n_part(self):
        # _calc_s_b_to_s_n_part
        for i in range(1, 7+1):
            ym, s_n, yd = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i))
            s_b = s_n - 1
            yd_s_n = copy.deepcopy(yd)
            _, bl_tuple = load_branching_law(s_n, yd_s_n, is_flag_true_if_not_s_n=False, is_return_tuple=True)
            bl_num, _, _, before_yd = bl_tuple
            _, yt_all_s_n = load_young_table(s_n, yd_s_n, is_flag_true_if_not_s_n=False)
            _, matrix_div = load_young_table_num(s_n, yd_s_n, is_flag_true_if_not_s_n=False)

            # _calc_s_b_to_s_n_part
            flag, matrix_in = _calc_s_b_to_s_n_part(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, yt_all_s_n)
            assert flag
            assert matrix_in == ym

    # def test_orthonormal(self):
    #     # TODO 测试正交归一性
    #     pass
    #
    # def test_trace(self):
    #     # TODO 测试矩阵的迹
    #     pass


# np的计算方式已经不再支持
# @pytest.mark.skip("pass")
# class TestYamanouchiMatrixNP(object):
#
#     def setup_class(self):
#         self.protector = DBProtector(cgc_db_name, extension_name=".test_yamanouchi_matrix_protected")
#         self.protector.protector_setup()
#
#         # 准备前文
#         flag, msg = create_young_diagrams(6)
#         assert flag
#         assert msg == 6
#         flag, msg = create_branching_laws(6)
#         assert flag
#         assert msg == 6
#         flag, data = create_young_tableaux(6)
#         assert flag
#         assert msg == 6
#
#         # (Sb, Sn)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
#         self.ym_1 = np.array([[-1/2, math.sqrt(3)/2], [math.sqrt(3)/2, 1/2]])
#         self.sn_1, self.yd_1, self.ix_1 = 3, [2, 1], (2, 3,)
#         self.ym_2 = np.array([[-1/3, math.sqrt(8)/3, 0], [math.sqrt(8)/3, 1/3, 0], [0, 0, 1]])
#         self.sn_2, self.yd_2, self.ix_2 = 4, [3, 1], (3, 4,)
#         self.ym_3 = np.array([[-1, 0, 0], [0, -1/3, math.sqrt(8)/3], [0, math.sqrt(8)/3, 1/3]])
#         self.sn_3, self.yd_3, self.ix_3 = 4, [2, 1, 1], (3, 4,)
#         self.ym_4 = np.array([[1, 0], [0, -1]])
#         self.sn_4, self.yd_4, self.ix_4 = 4, [2, 2], (3, 4,)
#         self.ym_5 = np.array([[-1/4, math.sqrt(15)/4, 0, 0], [math.sqrt(15)/4, 1/4, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
#         self.sn_5, self.yd_5, self.ix_5 = 5, [4, 1], (4, 5,)
#         self.ym_6 = np.array([[1, 0, 0, 0, 0], [0, -1/2, 0, math.sqrt(3)/2, 0], [0, 0, -1/2, 0, math.sqrt(3)/2],
#                               [0, math.sqrt(3)/2, 0, 1/2, 0], [0, 0, math.sqrt(3)/2, 0, 1/2]])
#         self.sn_6, self.yd_6, self.ix_6 = 5, [3, 2], (4, 5,)
#         self.ym_7 = np.array([[-1, 0, 0, 0, 0, 0], [0, -1/4, 0, math.sqrt(15)/4, 0, 0],
#                               [0, 0, -1/4, 0, math.sqrt(15)/4, 0], [0, math.sqrt(15)/4, 0, 1/4, 0, 0],
#                               [0, 0, math.sqrt(15)/4, 0, 1/4, 0], [0, 0, 0, 0, 0, 1]])
#         self.sn_7, self.yd_7, self.ix_7 = 5, [3, 1, 1], (4, 5,)
#
#         # (1, 2)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
#         self.ym_11 = np.array([[1, 0], [0, -1]])
#         self.sn_11, self.yd_11, self.ix_11 = 3, [2, 1], (1, 2,)
#         self.ym_12 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
#         self.sn_12, self.yd_12, self.ix_12 = 4, [3, 1], (1, 2,)
#         self.ym_13 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
#         self.sn_13, self.yd_13, self.ix_13 = 4, [2, 1, 1], (1, 2,)
#         self.ym_14 = np.array([[1, 0], [0, -1]])
#         self.sn_14, self.yd_14, self.ix_14 = 4, [2, 2], (1, 2,)
#         self.ym_15 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
#         self.sn_15, self.yd_15, self.ix_15 = 5, [4, 1], (1, 2,)
#         self.ym_16 = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, -1]])
#         self.sn_16, self.yd_16, self.ix_16 = 5, [3, 2], (1, 2,)
#         self.ym_17 = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0],
#                                [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]])
#         self.sn_17, self.yd_17, self.ix_17 = 5, [3, 1, 1], (1, 2,)
#
#         # (2, 3)  # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2
#         self.ym_21 = np.array([[-1/2, math.sqrt(3)/2], [math.sqrt(3)/2, 1/2]])
#         self.sn_21, self.yd_21, self.ix_21 = 3, [2, 1], (2, 3,)
#         self.ym_22 = np.array([[1, 0, 0], [0, -1/2, math.sqrt(3)/2], [0, math.sqrt(3)/2, 1/2]])
#         self.sn_22, self.yd_22, self.ix_22 = 4, [3, 1], (2, 3,)
#         self.ym_23 = np.array([[-1/2, math.sqrt(3)/2, 0], [math.sqrt(3)/2, 1/2, 0], [0, 0, -1]])
#         self.sn_23, self.yd_23, self.ix_23 = 4, [2, 1, 1], (2, 3,)
#         self.ym_24 = np.array([[-1/2, math.sqrt(3)/2], [math.sqrt(3)/2, 1/2]])
#         self.sn_24, self.yd_24, self.ix_24 = 4, [2, 2], (2, 3,)
#         self.ym_25 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1/2, math.sqrt(3)/2], [0, 0, math.sqrt(3)/2, 1/2]])
#         self.sn_25, self.yd_25, self.ix_25 = 5, [4, 1], (2, 3,)
#         self.ym_26 = np.array([[1, 0, 0, 0, 0], [0, -1/2, math.sqrt(3)/2, 0, 0], [0, math.sqrt(3)/2, 1/2, 0, 0],
#                                [0, 0, 0, -1/2, math.sqrt(3)/2], [0, 0, 0, math.sqrt(3)/2, 1/2]])
#         self.sn_26, self.yd_26, self.ix_26 = 5, [3, 2], (2, 3,)
#         self.ym_27 = np.array([[1, 0, 0, 0, 0, 0], [0, -1/2, math.sqrt(3)/2, 0, 0, 0],
#                                [0, math.sqrt(3)/2, 1/2, 0, 0, 0], [0, 0, 0, -1/2, math.sqrt(3)/2, 0],
#                                [0, 0, 0, math.sqrt(3)/2, 1/2, 0], [0, 0, 0, 0, 0, -1]])
#         self.sn_27, self.yd_27, self.ix_27 = 5, [3, 1, 1], (2, 3,)
#
#         # (1, n)  # 我之前算的
#         self.ym_31 = np.array([[-0.4999999999999999, -0.8660254037844386], [-0.8660254037844386, 0.4999999999999999]])
#         self.sn_31, self.yd_31, self.ix_31 = 3, [2, 1], (1, 3,)
#         self.ym_32 = np.array([[-0.33333333333333315, -0.4714045207910316, -0.8164965809277259],
#                                [-0.4714045207910316, 0.8333333333333331, -0.28867513459481287],
#                                [-0.8164965809277259, -0.28867513459481287, 0.4999999999999999]])
#         self.sn_32, self.yd_32, self.ix_32 = 4, [3, 1], (1, 4,)
#         self.ym_33 = np.array([[-0.4999999999999998, -0.2886751345948128, 0.8164965809277259],
#                                [-0.2886751345948128, -0.8333333333333334, -0.4714045207910316],
#                                [0.8164965809277259, -0.4714045207910316, 0.3333333333333332]])
#         self.sn_33, self.yd_33, self.ix_33 = 4, [2, 1, 1], (1, 4,)
#         self.ym_34 = np.array([[-0.4999999999999999, 0.8660254037844386], [0.8660254037844386, 0.4999999999999999]])
#         self.sn_34, self.yd_34, self.ix_34 = 4, [2, 2], (1, 4,)
#
#         # (2, n)  # 我之前算的
#         self.ym_41 = np.array([[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]])
#         self.sn_41, self.yd_41, self.ix_41 = 3, [2, 1], (2, 3,)
#         self.ym_42 = np.array([[-0.3333333333333333, -0.4714045207910317, 0.8164965809277259],
#                                [-0.4714045207910317, 0.8333333333333333, 0.28867513459481287],
#                                [0.8164965809277259, 0.28867513459481287, 0.5]])
#         self.sn_42, self.yd_42, self.ix_42 = 4, [3, 1], (2, 4,)
#         self.ym_43 = np.array([[-0.5, 0.28867513459481287, -0.8164965809277259],
#                                [0.28867513459481287, -0.8333333333333334, -0.4714045207910317],
#                                [-0.8164965809277259, -0.4714045207910317, 0.3333333333333333]])
#         self.sn_43, self.yd_43, self.ix_43 = 4, [2, 1, 1], (2, 4,)
#         self.ym_44 = np.array([[-0.5, -0.8660254037844386], [-0.8660254037844386, 0.5]])
#         self.sn_44, self.yd_44, self.ix_44 = 4, [2, 2], (2, 4,)
#         pass
#
#     def teardown_class(self):
#         self.protector.protector_teardown()
#         pass
#
#     # start with 0xx tests need test by order
#
#     def test_001_create_yamanouchi_matrix_sn_2(self):
#         # check with no db
#         flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
#         assert flag
#         assert finish_s_n == 1
#
#         # create
#         flag, test_s_n = create_yamanouchi_matrix(2, is_sympy=False)
#         assert flag
#         assert test_s_n == 2
#
#         # check answers
#         all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
#         for i in all_ij_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#         all_in_num_list = list(range(1, 7 + 1)) + list(range(31, 34 + 1)) + list(range(41, 44 + 1))
#         for i in all_in_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#
#         # check finish s_n
#         flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
#         assert flag
#         assert finish_s_n == 2
#
#         # history_times
#         _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
#         flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
#         assert flag
#         assert data.get("data") == []
#         assert isinstance(data.get("flags"), dict)
#         assert isinstance(data.get("flags").get("history_times"), dict)
#         for i in [2]:
#             assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
#             assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
#         flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
#         assert isinstance(data_txt, str)
#         data = eval(data_txt)
#         assert isinstance(data, dict)
#         for i in [2]:
#             assert isinstance(data.get("history_times").get("S{}".format(i)), int)
#
#     def test_002_create_yamanouchi_matrix_sn_3_to_4(self):
#         # create
#         flag, test_s_n = create_yamanouchi_matrix(4, is_sympy=False)
#         assert flag
#         assert test_s_n == 4
#
#         # check answers
#         all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
#         for i in all_ij_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#         all_in_num_list = list(range(1, 7 + 1)) + list(range(31, 34 + 1)) + list(range(41, 44 + 1))
#         for i in all_in_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#
#         # check finish s_n
#         flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
#         assert flag
#         assert finish_s_n == 4
#
#         # history_times
#         _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
#         flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
#         assert flag
#         assert data.get("data") == []
#         assert isinstance(data.get("flags"), dict)
#         assert isinstance(data.get("flags").get("history_times"), dict)
#         for i in range(2, test_s_n + 1):
#             assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
#             assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
#         flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
#         assert isinstance(data_txt, str)
#         data = eval(data_txt)
#         assert isinstance(data, dict)
#         for i in range(2, test_s_n + 1):
#             assert isinstance(data.get("history_times").get("S{}".format(i)), int)
#
#     def test_002_create_yamanouchi_matrix_sn_5_to_6(self):
#         # create
#         flag, test_s_n = create_yamanouchi_matrix(6, is_sympy=False)
#         assert flag
#         assert test_s_n == 6
#
#         # check answers
#         all_ij_num_list = list(range(1, 7 + 1)) + list(range(11, 17 + 1)) + list(range(21, 27 + 1))
#         for i in all_ij_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#         all_in_num_list = list(range(1, 7 + 1)) + list(range(31, 34 + 1)) + list(range(41, 44 + 1))
#         for i in all_in_num_list:
#             ym, s_i = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i))
#             ix, yd = eval("self.ix_{}".format(i)), eval("self.yd_{}".format(i))
#             if s_i <= test_s_n:
#                 _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
#                 flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
#                 assert flag
#                 assert isinstance(data.get("create_time"), str)
#                 assert np.allclose(data.get("data"), ym, atol=0.0000000000001), "s_i={}, yd={}, ix={}, i={}".format(
#                     s_i, yd, ix, i)
#                 assert data.get("flags").get("speed_time") >= 0
#
#         # check finish s_n
#         flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
#         assert flag
#         assert finish_s_n == 6
#
#         # history_times
#         _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
#         flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
#         assert flag
#         assert data.get("data") == []
#         assert isinstance(data.get("flags"), dict)
#         assert isinstance(data.get("flags").get("history_times"), dict)
#         for i in range(2, test_s_n + 1):
#             assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
#             assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
#         flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
#         assert isinstance(data_txt, str)
#         data = eval(data_txt)
#         assert isinstance(data, dict)
#         for i in range(2, test_s_n + 1):
#             assert isinstance(data.get("history_times").get("S{}".format(i)), int)
#
#     def test__calc_s_b_to_s_n_part(self):
#         # _calc_s_b_to_s_n_part
#         for i in range(1, 7+1):
#             ym, s_n, yd = eval("self.ym_{}".format(i)), eval("self.sn_{}".format(i)), eval("self.yd_{}".format(i))
#             s_b = s_n - 1
#             yd_s_n = copy.deepcopy(yd)
#             _, bl_tuple = load_branching_law(s_n, yd_s_n, is_flag_true_if_not_s_n=False, is_return_tuple=True)
#             bl_num, _, _, before_yd = bl_tuple
#             _, yt_all_s_n = load_young_table(s_n, yd_s_n, is_flag_true_if_not_s_n=False)
#             _, matrix_div = load_young_table_num(s_n, yd_s_n, is_flag_true_if_not_s_n=False)
#
#             # _calc_s_b_to_s_n_part
#             flag, matrix_in = _calc_s_b_to_s_n_part(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, yt_all_s_n,
#                                                     is_sympy=False)
#             assert flag
#             assert np.allclose(ym, matrix_in, atol=0.0000000000001)
#
#     def test_orthonormal(self):
#         # TODO 测试正交归一性
#         pass
#
#     def test_trace(self):
#         # TODO 测试矩阵的迹
#         pass
