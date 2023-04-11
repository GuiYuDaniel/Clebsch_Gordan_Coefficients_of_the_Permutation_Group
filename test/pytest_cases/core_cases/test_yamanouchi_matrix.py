# -*- coding:utf8 -*-
"""
测试
core/yamanouchi_matrix.py
下所有功能是否正确执行
"""


import copy
import sympy as sp
from sympy import Matrix, sqrt
from sympy import Rational as Ra
import pytest
from conf.cgc_config import cgc_db_name
from core.cgc_utils.cgc_db_typing import YamanouchiMatrixInfo
from core.branching_laws import create_branching_laws
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from core.young_tableaux import create_young_tableaux
from core.yamanouchi_matrix import create_yamanouchi_matrix
from core.yamanouchi_matrix import get_yamanouchi_matrix_finish_s_n
from core.cgc_utils.cgc_local_db import get_yamanouchi_matrix_file_name, get_yamanouchi_matrix_finish_s_n_name
from core.yamanouchi_matrix import _calc_s_b_to_s_n_part, load_yamanouchi_matrix
from core.branching_laws import load_branching_law
from core.young_tableaux import load_young_table, load_young_table_num
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class YamanouchiMatrixData(object):

    def __init__(self):
        # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2(全表都检查)
        # S3
        self.sn_31, self.yd_31, self.ym_31_dict = 3, [2, 1], {
            (1, 2): Matrix([[1, 0], [0, -1]]),
            (2, 3): Matrix([[Ra(-1) / 2, sqrt(3) / 2], [sqrt(3) / 2, Ra(1) / 2]])
        }

        # S4
        self.sn_41, self.yd_41, self.ym_41_dict = 4, [3, 1], {
            (1, 2): Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),
            (2, 3): Matrix([[1, 0, 0], [0, Ra(-1) / 2, sqrt(3) / 2], [0, sqrt(3) / 2, Ra(1) / 2]]),
            (3, 4): Matrix([[Ra(-1) / 3, sqrt(8) / 3, 0], [sqrt(8) / 3, Ra(1) / 3, 0], [0, 0, 1]])
        }
        self.sn_42, self.yd_42, self.ym_42_dict = 4, [2, 1, 1], {
            (1, 2): Matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),
            (2, 3): Matrix([[Ra(-1) / 2, sqrt(3) / 2, 0], [sqrt(3) / 2, Ra(1) / 2, 0], [0, 0, -1]]),
            (3, 4): Matrix([[-1, 0, 0], [0, Ra(-1) / 3, sqrt(8) / 3], [0, sqrt(8) / 3, Ra(1) / 3]])
        }
        self.sn_43, self.yd_43, self.ym_43_dict = 4, [2, 2], {
            (1, 2): Matrix([[1, 0], [0, -1]]),
            (2, 3): Matrix([[Ra(-1) / 2, sqrt(3) / 2], [sqrt(3) / 2, Ra(1) / 2]]),
            (3, 4): Matrix([[1, 0], [0, -1]])
        }

        # S5
        self.sn_51, self.yd_51, self.ym_51_dict = 5, [4, 1], {
            (1, 2): Matrix([[1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, -1]]),
            (2, 3): Matrix([[1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, Ra(-1) / 2, sqrt(3) / 2],
                            [0, 0, sqrt(3) / 2, Ra(1) / 2]]),
            (3, 4): Matrix([[1, 0, 0, 0],
                            [0, Ra(-1) / 3, sqrt(8) / 3, 0],
                            [0, sqrt(8) / 3, Ra(1) / 3, 0],
                            [0, 0, 0, 1]]),
            (4, 5): Matrix([[Ra(-1) / 4, sqrt(15) / 4, 0, 0],
                            [sqrt(15) / 4, Ra(1) / 4, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]])
        }
        self.sn_52, self.yd_52, self.ym_52_dict = 5, [3, 2], {
            (1, 2): Matrix([[1, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0],
                            [0, 0, -1, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, -1]]),
            (2, 3): Matrix([[1, 0, 0, 0, 0],
                            [0, Ra(-1) / 2, sqrt(3) / 2, 0, 0],
                            [0, sqrt(3) / 2, Ra(1) / 2, 0, 0],
                            [0, 0, 0, Ra(-1) / 2, sqrt(3) / 2],
                            [0, 0, 0, sqrt(3) / 2, Ra(1) / 2]]),
            (3, 4): Matrix([[Ra(-1) / 3, sqrt(8) / 3, 0, 0, 0],
                            [sqrt(8) / 3, Ra(1) / 3, 0, 0, 0],
                            [0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, -1]]),
            (4, 5): Matrix([[1, 0, 0, 0, 0],
                            [0, Ra(-1) / 2, 0, sqrt(3) / 2, 0],
                            [0, 0, Ra(-1) / 2, 0, sqrt(3) / 2],
                            [0, sqrt(3) / 2, 0, Ra(1) / 2, 0],
                            [0, 0, sqrt(3) / 2, 0, Ra(1) / 2]])
        }
        self.sn_53, self.yd_53, self.ym_53_dict = 5, [3, 1, 1], {
            (1, 2): Matrix([[1, 0, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0, 0],
                            [0, 0, -1, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0],
                            [0, 0, 0, 0, -1, 0],
                            [0, 0, 0, 0, 0, -1]]),
            (2, 3): Matrix([[1, 0, 0, 0, 0, 0],
                            [0, Ra(-1) / 2, sqrt(3) / 2, 0, 0, 0],
                            [0, sqrt(3) / 2, Ra(1) / 2, 0, 0, 0],
                            [0, 0, 0, Ra(-1) / 2, sqrt(3) / 2, 0],
                            [0, 0, 0, sqrt(3) / 2, Ra(1) / 2, 0],
                            [0, 0, 0, 0, 0, -1]]),
            (3, 4): Matrix([[Ra(-1) / 3, sqrt(8) / 3, 0, 0, 0, 0],
                            [sqrt(8) / 3, Ra(1) / 3, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0],
                            [0, 0, 0, -1, 0, 0],
                            [0, 0, 0, 0, Ra(-1) / 3, sqrt(8) / 3],
                            [0, 0, 0, 0, sqrt(8) / 3, Ra(1) / 3]]),
            (4, 5): Matrix([[-1, 0, 0, 0, 0, 0],
                            [0, Ra(-1) / 4, 0, sqrt(15) / 4, 0, 0],
                            [0, 0, -1 / 4, 0, sqrt(15) / 4, 0],
                            [0, sqrt(15) / 4, 0, Ra(1) / 4, 0, 0],
                            [0, 0, sqrt(15) / 4, 0, Ra(1) / 4, 0],
                            [0, 0, 0, 0, 0, 1]])
        }
        self.sn_54, self.yd_54, self.ym_54_dict = 5, [2, 2, 1], {
            (1, 2): Matrix([[1, 0, 0, 0, 0],
                            [0, -1, 0, 0, 0],
                            [0, 0, 1, 0, 0],
                            [0, 0, 0, -1, 0],
                            [0, 0, 0, 0, -1]]),
            (2, 3): Matrix([[Ra(-1) / 2, sqrt(3) / 2, 0, 0, 0],
                            [sqrt(3) / 2, Ra(1) / 2, 0, 0, 0],
                            [0, 0, Ra(-1) / 2, sqrt(3) / 2, 0],
                            [0, 0, sqrt(3) / 2, Ra(1) / 2, 0],
                            [0, 0, 0, 0, -1]]),
            (3, 4): Matrix([[1, 0, 0, 0, 0],
                            [0, -1, 0, 0, 0],
                            [0, 0, -1, 0, 0],
                            [0, 0, 0, Ra(-1) / 3, sqrt(8) / 3],
                            [0, 0, 0, sqrt(8) / 3, Ra(1) / 3]]),
            (4, 5): Matrix([[Ra(-1) / 2, 0, sqrt(3) / 2, 0, 0],
                            [0, Ra(-1) / 2, 0, sqrt(3) / 2, 0],
                            [sqrt(3) / 2, 0, Ra(1) / 2, 0, 0],
                            [0, sqrt(3) / 2, 0, Ra(1) / 2, 0],
                            [0, 0, 0, 0, -1]])
        }
        self.sn_55, self.yd_55, self.ym_55_dict = 5, [2, 1, 1, 1], {
            (1, 2): Matrix([[1, 0, 0, 0],
                            [0, -1, 0, 0],
                            [0, 0, -1, 0],
                            [0, 0, 0, -1]]),
            (2, 3): Matrix([[Ra(-1) / 2, sqrt(3) / 2, 0, 0],
                            [sqrt(3) / 2, Ra(1) / 2, 0, 0],
                            [0, 0, -1, 0],
                            [0, 0, 0, -1]]),
            (3, 4): Matrix([[-1, 0, 0, 0],
                            [0, Ra(-1) / 3, sqrt(8) / 3, 0],
                            [0, sqrt(8) / 3, Ra(1) / 3, 0],
                            [0, 0, 0, -1]]),
            (4, 5): Matrix([[-1, 0, 0, 0],
                            [0, -1, 0, 0],
                            [0, 0, Ra(-1) / 4, sqrt(15) / 4],
                            [0, 0, sqrt(15) / 4, Ra(1) / 4]])
        }

        self.ym_num_list = [31, 41, 42, 43, 51, 52, 53, 54, 55]


# @pytest.mark.skip("pass")
class TestYamanouchiMatrixSP(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_db_name, extension_name=".test_yamanouchi_matrix_protected")
        self.protector.protector_setup()

        # 准备前文
        self.test_sn = 7
        flag, msg = create_young_diagrams(self.test_sn)
        assert flag
        assert msg == self.test_sn
        flag, msg = create_branching_laws(self.test_sn)
        assert flag
        assert msg == self.test_sn
        flag, data = create_young_tableaux(self.test_sn)
        assert flag
        assert msg == self.test_sn

        # 全部来自《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节 表4.4-2(全表都检查)
        self.data = YamanouchiMatrixData()
        self.ym_num_list = self.data.ym_num_list
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
        start_sn = 2
        end_sn = 2
        flag, test_s_n = create_yamanouchi_matrix(end_sn)
        assert flag
        assert test_s_n == end_sn

        # check answers
        for num in self.ym_num_list:
            s_i, yd = eval("self.data.sn_{}".format(num)), eval("self.data.yd_{}".format(num))
            ym_dict = eval("self.data.ym_{}_dict".format(num))
            if s_i > end_sn or s_i < start_sn:
                continue
            for ix, ym in ym_dict.items():
                if ix[1] - ix[0] == 1:  # check (i,j)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, ix={}, rst={}".format(num, ix, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0
                if ix[1] == s_i:  # check (i,n)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, ix={}, rst={}".format(num, ix, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == end_sn

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_002_create_yamanouchi_matrix_sn_3_to_4(self):
        # create
        start_sn = 3
        end_sn = 4
        flag, test_s_n = create_yamanouchi_matrix(end_sn)
        assert flag
        assert test_s_n == end_sn

        # check answers
        for num in self.ym_num_list:
            s_i, yd = eval("self.data.sn_{}".format(num)), eval("self.data.yd_{}".format(num))
            ym_dict = eval("self.data.ym_{}_dict".format(num))
            if s_i > end_sn or s_i < start_sn:
                continue
            for ix, ym in ym_dict.items():
                if ix[1] - ix[0] == 1:  # check (i,j)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, rst={}".format(num, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0
                if ix[1] == s_i:  # check (i,n)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, rst={}".format(num, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == end_sn

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test_003_create_yamanouchi_matrix_sn_5_to_7(self):
        # create
        start_sn = 5
        end_sn = 7
        flag, test_s_n = create_yamanouchi_matrix(end_sn)
        assert flag
        assert test_s_n == end_sn

        # check answers
        for num in self.ym_num_list:
            s_i, yd = eval("self.data.sn_{}".format(num)), eval("self.data.yd_{}".format(num))
            ym_dict = eval("self.data.ym_{}_dict".format(num))
            if s_i > end_sn or s_i < start_sn:
                continue
            for ix, ym in ym_dict.items():
                if ix[1] - ix[0] == 1:  # check (i,j)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="ij")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, ix={}, rst={}".format(num, ix, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0
                if ix[1] == s_i:  # check (i,n)
                    _, file_name = get_yamanouchi_matrix_file_name(s_i, yd, ix, mode="in")
                    flag, data = YamanouchiMatrixInfo(s_i).query_by_file_name(file_name)
                    assert flag
                    assert isinstance(data.get("create_time"), str)
                    assert data.get("data") == ym, "num={}, ix={}, rst={}".format(num, ix, data.get("data"))
                    assert data.get("flags").get("speed_time") >= 0

        # check finish s_n
        flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
        assert flag
        assert finish_s_n == end_sn

        # history_times
        _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
        flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == []
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 20
        flag, data_txt = YamanouchiMatrixInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(start_sn, end_sn + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    def test__calc_s_b_to_s_n_part(self):
        # _calc_s_b_to_s_n_part
        for num in self.ym_num_list:
            s_i, yd = eval("self.data.sn_{}".format(num)), eval("self.data.yd_{}".format(num))
            s_b = s_i - 1
            ym_dict = eval("self.data.ym_{}_dict".format(num))
            for ix, ym in ym_dict.items():
                if ix[0] == s_b and ix[1] == s_i:
                    yd_s_n = copy.deepcopy(yd)
                    _, bl_tuple = load_branching_law(s_i, yd_s_n, is_flag_true_if_not_s_n=False, is_return_tuple=True)
                    bl_num, _, _, before_yd = bl_tuple
                    _, yt_all_s_n = load_young_table(s_i, yd_s_n, is_flag_true_if_not_s_n=False)
                    _, matrix_div = load_young_table_num(s_i, yd_s_n, is_flag_true_if_not_s_n=False)

                    # _calc_s_b_to_s_n_part
                    flag, matrix_in = _calc_s_b_to_s_n_part(matrix_div, s_b, s_i, yd_s_n, bl_num, before_yd, yt_all_s_n)
                    assert flag
                    assert matrix_in == ym


# @pytest.mark.skip("pass")
class TestResultsWithoutCalc(object):

    def setup_class(self):
        # 它不改变数据库，所以不需要保护

        # 检查范围
        self.start_sn = 2
        self.end_sn = 7
        Matrix().is_symmetric()

    def teardown_class(self):
        pass

    def test_orthonormal_and_trace(self):
        """利用施密特正交化测正交归一性，测试矩阵的迹"""
        # 取数据
        for s_i in range(self.start_sn, self.end_sn + 1):
            ix_dict = {"ij": [(i + 1, i + 2) for i in range(s_i - 1)],
                       "in": [(i + 1, s_i) for i in range(s_i - 1)]}
            _, yds = load_young_diagrams(s_i)
            for yd in yds:
                first_yd_trace = None  # 每个yd的迹是可以计算的，但后续用不到，就没写。这里，只需对比同yd，trace相等即可
                for ix_key, ix_list in ix_dict.items():
                    for ix in ix_list:
                        _, ym = load_yamanouchi_matrix(s_i, yd, ix, mode=ix_key)
                        # 检查对称性（实厄米，就是对称）
                        assert ym.is_symmetric()
                        # 检查正交归一性，因为对称，行/列查一个就行
                        ym_by_row = [ym.row(i) for i in range(ym.rows)]
                        assert ym_by_row == sp.GramSchmidt(ym_by_row, True)
                        # 正交归一的对称矩阵独有特性!!!
                        assert sp.eye(ym.rows) == ym * ym
                        # 矩阵的迹
                        if first_yd_trace is None:
                            first_yd_trace = ym.trace()
                        else:
                            assert first_yd_trace == ym.trace()
