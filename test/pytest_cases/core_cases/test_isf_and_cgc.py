# -*- coding:utf8 -*-
"""
测试
core/isf_and_cgc.py
下所有功能是否正确执行
"""


import os
import time
import math
import numpy as np
import sympy as sp
from sympy import Rational as Ra
from sympy import sqrt
import pytest
from itertools import chain, product
from conf.cgc_config import cgc_rst_folder, isf_0_error_value
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from core.branching_laws import create_branching_laws
from core.young_tableaux import create_young_tableaux
from core.yamanouchi_matrix import create_yamanouchi_matrix
from core.characters_and_gi import create_characters_and_gi
from core.cg_series import create_cg_series
from core.eigenvalues import create_eigenvalues, load_eigenvalues
from core.cgc_utils.cgc_db_typing import ISFInfo, CGCInfo
from core.cgc_utils.cgc_local_db import get_isf_file_name, get_isf_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_cgc_file_name, get_cgc_finish_s_n_name
from core.isf_and_cgc import create_isf_and_cgc, load_isf, load_cgc, load_cgc_with_m1
from core.isf_and_cgc import get_isf_finish_s_n, get_cgc_finish_s_n
from core.isf_and_cgc import ΣMDataHelper, DataHelper, CalcHelper, ISFHelper, CGCHelper
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class Helper(object):
    """这里放一些非直接测试的辅助函数"""

    def __init__(self, s_t_max=5):
        # 这里拿yd对应的eigenvalues_list是为了判断(σ', μ')的Yamanouchi顺序，在S5之前，它都是和Yamanouchi序一致的
        self.all_yd_dict = {}
        self.all_ev_dict = {}
        for s_t in range(1, s_t_max + 1):
            # yd
            _, yd_list = load_young_diagrams(s_t, is_flag_true_if_not_s_n=False)
            self.all_yd_dict[s_t] = yd_list

            # eigenvalue
            _, eigenvalue_list = load_eigenvalues(s_t, is_flag_true_if_not_s_n=False)
            self.all_ev_dict[s_t] = eigenvalue_list

    def calc_symmetry_cgc(self, s_n, σ, μ, ν, β, m, cgc_square_dict):
        """CGC_{μm_μ,σm_σ,νβm} = ϵ(σμνβ) * CGC_{σm_σ,μm_μ,νβm}

        return和Data中self.cgc_x格式一致"""
        ϵ = self.calc_ϵ(s_n, σ, μ, ν, β)
        symmetry_cgc_square_dict = {(k[1], k[0]): ϵ * v for k, v in cgc_square_dict.items() if k != "N"}
        symmetry_cgc_square_dict["N"] = cgc_square_dict["N"]
        return s_n, μ, σ, ν, β, m, symmetry_cgc_square_dict

    def _calc_symmetry_isf_rows(self, s_n, isf_rows):
        all_yd_list = self.all_yd_dict[s_n - 1]
        all_ev_list = self.all_ev_dict[s_n - 1]
        symmetry_isf_rows = []
        μ_st_ev_set = set(all_ev_list[all_yd_list.index(i[1])] for i in isf_rows)
        μ_st_ev_list = list(μ_st_ev_set)
        for μ_st_ev in sorted(μ_st_ev_list, reverse=True):
            μ_st = all_yd_list[all_ev_list.index(μ_st_ev)]
            σ_st_ev_list_part = [all_ev_list[all_yd_list.index(i[0])] for i in isf_rows if i[1] == μ_st]
            descending_σ_st_ev_list_part = sorted(σ_st_ev_list_part, reverse=True)
            symmetry_isf_rows_part = [(μ_st, all_yd_list[all_ev_list.index(σ_st_ev)])
                                      for σ_st_ev in descending_σ_st_ev_list_part]
            symmetry_isf_rows += symmetry_isf_rows_part
        return symmetry_isf_rows

    @staticmethod
    def _calc_symmetry_isf_square_matrix_tmp(isf_rows, symmetry_isf_rows, isf_square_matrix):
        isf_anti_rows = [(i[1], i[0]) for i in isf_rows]  # S6之前，σ', μ'都可用  # 它是简单调换顺序，用来对比看哪些行移动了的
        # 对比symmetry_isf_rows 和 isf_anti_rows，得出哪些行应该变动
        symmetry_isf_tmp = sp.zeros(len(isf_rows))  # 还可以用交换来优化
        for i, row in enumerate(symmetry_isf_rows):
            symmetry_isf_tmp[i, :] = isf_square_matrix[isf_anti_rows.index(row), :]
        return symmetry_isf_tmp

    def calc_symmetry_isf_square_dict(self, s_n, isf_square_dict, σ, μ, ν_st):
        isf_rows = isf_square_dict["rows"]  # rows需要按照σ，μ调换顺序，不仅仅左右换位置，还要考虑到Yamanouchi序
        symmetry_isf_cols = isf_square_dict["cols"]  # cols是不变的
        symmetry_isf_rows = self._calc_symmetry_isf_rows(s_n, isf_rows)
        ϵ_list = []
        for col in symmetry_isf_cols:
            ν, β = (col, None) if isinstance(col, list) else (col[0], col[1])
            ϵ = self.calc_ϵ(s_n, σ, μ, ν, β)
            ϵ_list.append(ϵ)
        ϵ_st_list = []
        for row in symmetry_isf_rows:
            μ_st, σ_st, β_st = (row[0], row[1], None) if len(row) == 2 else (row[0], row[1], row[2])
            ϵ_st = self.calc_ϵ(s_n-1, σ_st, μ_st, ν_st, β_st)
            ϵ_st_list.append(ϵ_st)
        isf_square_matrix = isf_square_dict["isf"]
        symmetry_isf_tmp = self._calc_symmetry_isf_square_matrix_tmp(isf_rows, symmetry_isf_rows, isf_square_matrix)
        div = len(isf_rows)
        symmetry_isf_square = sp.zeros(div)
        for (row_i, ϵ_st), (col_i, ϵ) in product(zip(range(div), ϵ_st_list), zip(range(div), ϵ_list)):
            symmetry_isf_square[row_i, col_i] = ϵ * ϵ_st * symmetry_isf_tmp[row_i, col_i]
        return {"rows": symmetry_isf_rows,
                "cols": symmetry_isf_cols,
                "isf": symmetry_isf_square}

    def calc_symmetry_isf(self, s_n, σ, μ, ν_st, isf_square_dict):
        symmetry_isf_dict = self.calc_symmetry_isf_square_dict(s_n, isf_square_dict, σ, μ, ν_st)
        return s_n, μ, σ, ν_st, symmetry_isf_dict

    @staticmethod
    def get_min_μ_σ_from_dict(cgc_square_dict):
        m_μ_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple))
        min_m_μ = min(m_μ_set)
        m_σ_set_of_min_m_μ = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[1] == min_m_μ)
        min_m_σ = min(m_σ_set_of_min_m_μ)
        min_key = (min_m_σ, min_m_μ,)
        return min_key, cgc_square_dict[min_key]

    def calc_ϵ(self, s_n, σ, μ, ν, β):
        """ϵ(σμνβ) = sign(CGC_{σm_σ,μm_μ,νβm_1} | (m_μ,m_σ)=min)"""
        flag, cgc_square_dict = load_cgc_with_m1(s_n, σ, μ, ν, β)
        assert flag
        assert cgc_square_dict
        _, min_value = self.get_min_μ_σ_from_dict(cgc_square_dict)
        return sp.sign(min_value)

    @staticmethod
    def is_isf_square_orthogonalization(isf_square_matrix):
        """看一个平方方阵，所有行/列是否正交归一"""
        assert isf_square_matrix.shape[0] == isf_square_matrix.shape[1]
        div = isf_square_matrix.shape[0]
        for i in range(div):
            matrix_i_by_row = isf_square_matrix[i, :]
            matrix_i_by_col = isf_square_matrix[:, i]
            for j in range(div):
                matrix_j_by_row = isf_square_matrix[j, :]
                matrix_j_by_col = isf_square_matrix[:, j]

                sum_row_row = sum(sp.sign(matrix_i_by_row[k])
                                  * sp.sign(matrix_j_by_row[k])
                                  * sqrt(abs(matrix_i_by_row[k]))
                                  * sqrt(abs(matrix_j_by_row[k])) for k in range(div))
                sum_col_col = sum(sp.sign(matrix_i_by_col[k])
                                  * sp.sign(matrix_j_by_col[k])
                                  * sqrt(abs(matrix_i_by_col[k]))
                                  * sqrt(abs(matrix_j_by_col[k])) for k in range(div))
                if i == j:
                    assert sum_row_row == 1, "sum_row_row={} not eq 1 for i={},j={}".format(sum_row_row, i, j)
                    assert sum_col_col == 1, "sum_col_col={} not eq 1 for i={},j={}".format(sum_col_col, i, j)
                else:
                    assert sum_row_row == 0, "sum_row_row={} not eq 0 for i={},j={}".format(sum_row_row, i, j)
                    assert sum_col_col == 0, "sum_col_col={} not eq 0 for i={},j={}".format(sum_col_col, i, j)


class Data(object):
    """ISF 和 CGC 等数据，独立出来易读性好，让test类继承它"""

    def __init__(self):
        # isf data
        # 格式 Sn, σ, μ, ν_st, ISF_square_dict
        isf_square_dict = {"rows": [([1], [1])], "cols": [[2]], "isf": sp.Matrix([[1]])}
        self.isf_1 = (2, [2], [2], [1], isf_square_dict)

        # 数据来源见《群表示论的新途径》陈金全（上海科学技术出版社1984）表4.19
        isf_square_dict = {"rows": [([2], [2]), ([1, 1], [1, 1])],
                           "cols": [[3], [2, 1]],
                           "isf": sp.Matrix([[Ra(1) / 2, Ra(1) / 2], [Ra(1) / 2, -Ra(1) / 2]])}
        self.isf_2 = (3, [2, 1], [2, 1], [2], isf_square_dict)

        isf_square_dict = {"rows": [([2], [1, 1]), ([1, 1], [2])],
                           "cols": [[2, 1], [1, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 2, Ra(1) / 2], [-Ra(1) / 2, -Ra(1) / 2]])}
        self.isf_3 = (3, [2, 1], [2, 1], [1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3], [3]), ([2, 1], [2, 1])],
                           "cols": [[4], [3, 1]],
                           "isf": sp.Matrix([[Ra(1) / 3, Ra(2) / 3], [Ra(2) / 3, -Ra(1) / 3]])}
        self.isf_4 = (4, [3, 1], [3, 1], [3], isf_square_dict)

        isf_square_dict = {"rows": [([3], [2, 1]), ([2, 1], [3]), ([2, 1], [2, 1])],
                           "cols": [[3, 1], [2, 2], [2, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 6, Ra(1) / 3, Ra(1) / 2],
                                             [-Ra(1) / 6, Ra(1) / 3, -Ra(1) / 2],
                                             [Ra(2) / 3, Ra(1) / 3, 0]])}
        self.isf_5 = (4, [3, 1], [3, 1], [2, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1], [2, 1])], "cols": [[2, 1, 1]], "isf": sp.Matrix([[1]])}
        self.isf_6 = (4, [3, 1], [3, 1], [1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1], [2, 1])], "cols": [[2, 1, 1]], "isf": sp.Matrix([[-1]])}
        self.isf_7 = (4, [3, 1], [2, 2], [1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3], [2, 1]), ([2, 1], [2, 1]), ([2, 1], [1, 1, 1])],
                           "cols": [[3, 1], [2, 2], [2, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 2, Ra(1) / 3, Ra(1) / 6],
                                             [0, -Ra(1) / 3, Ra(2) / 3],
                                             [Ra(1) / 2, Ra(1) / 3, Ra(1) / 6]])}
        self.isf_8 = (4, [3, 1], [2, 1, 1], [2, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1], [2, 1]), ([2, 1], [1, 1, 1])],
                           "cols": [[3, 1], [2, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 2, Ra(1) / 2], [-Ra(1) / 2, -Ra(1) / 2]])}
        self.isf_9 = (4, [2, 2], [2, 1, 1], [2, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1], [2, 1]), ([1, 1, 1], [1, 1, 1])],
                           "cols": [[4], [3, 1]],
                           "isf": sp.Matrix([[Ra(2) / 3, Ra(1) / 3], [Ra(1) / 3, -Ra(2) / 3]])}
        self.isf_10 = (4, [2, 1, 1], [2, 1, 1], [3], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1], [2, 1]), ([2, 1], [1, 1, 1]), ([1, 1, 1], [2, 1])],
                           "cols": [[3, 1], [2, 2], [2, 1, 1]],
                           "isf": sp.Matrix([[Ra(2) / 3, Ra(1) / 3, 0],
                                             [-Ra(1) / 6, Ra(1) / 3, Ra(1) / 2],
                                             [-Ra(1) / 6, Ra(1) / 3, -Ra(1) / 2]])}
        self.isf_11 = (4, [2, 1, 1], [2, 1, 1], [2, 1], isf_square_dict)  # 这个也可以用来检查first_no_0

        isf_square_dict = {"rows": [([4], [3, 1]), ([3, 1], [4]), ([3, 1], [3, 1])],
                           "cols": [[4, 1], [3, 2], [3, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 12, Ra(5) / 12, Ra(1) / 2],
                                             [-Ra(1) / 12, Ra(5) / 12, -Ra(1) / 2],
                                             [Ra(5) / 6, Ra(1) / 6, 0]])}
        self.isf_12 = (5, [4, 1], [4, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([4], [2, 2]), ([3, 1], [3, 1])],
                           "cols": [[3, 2], [2, 2, 1]],
                           "isf": sp.Matrix([[-Ra(3) / 8, Ra(5) / 8], [-Ra(5) / 8, -Ra(3) / 8]])}
        self.isf_13 = (5, [4, 1], [3, 2], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1])],
                           "cols": [[3, 2], [2, 2, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 16, Ra(15) / 16], [Ra(15) / 16, Ra(1) / 16]])}
        self.isf_14 = (5, [4, 1], [3, 1, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1])],
                           "cols": [[2, 1, 1, 1]],
                           "isf": sp.Matrix([[-1]])}
        self.isf_15 = (5, [4, 1], [2, 2, 1], [1, 1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1])], "cols": [[3, 1, 1]], "isf": sp.Matrix([[1]])}
        self.isf_16 = (5, [4, 1], [2, 1, 1, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 2]), ([2, 2], [3, 1])],
                           "cols": [[4, 1], [3, 2], [3, 1, 1]],
                           "isf": sp.Matrix([[Ra(1) / 3, Ra(2) / 3, 0],
                                             [-Ra(1) / 3, Ra(1) / 6, Ra(1) / 2],
                                             [-Ra(1) / 3, Ra(1) / 6, -Ra(1) / 2]])}
        self.isf_17 = (5, [3, 2], [3, 2], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 2]), ([2, 2], [3, 1])],
                           "cols": [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]],
                           "isf": sp.Matrix([[Ra(2) / 5, 0, Ra(3) / 5],
                                             [Ra(3) / 10, Ra(1) / 2, -Ra(1) / 5],
                                             [-Ra(3) / 10, Ra(1) / 2, Ra(1) / 5]])}
        self.isf_18 = (5, [3, 2], [3, 2], [2, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 2], [2, 2])],
                           "cols": [[2, 1, 1, 1]],
                           "isf": sp.Matrix([[1]])}
        self.isf_19 = (5, [3, 2], [3, 2], [1, 1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1])],
                           "cols": [[4, 1]],
                           "isf": sp.Matrix([[1]])}
        self.isf_20 = (5, [3, 2], [3, 1, 1], [4], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1]), ([2, 2], [3, 1]), ([2, 2], [2, 1, 1])],
                           "cols": [[4, 1], [3, 2], ([3, 1, 1], 1), ([3, 1, 1], 2)],
                           "isf": sp.Matrix([[-Ra(3) / 10, 0, Ra(3) / 5, Ra(1) / 10],
                                             [-Ra(1) / 6, Ra(1) / 3, 0, -Ra(1) / 2],
                                             [-Ra(1) / 30, Ra(5) / 12, -Ra(3) / 20, Ra(2) / 5],
                                             [Ra(1) / 2, Ra(1) / 4, Ra(1) / 4, 0]])}
        self.isf_21 = (5, [3, 2], [3, 1, 1], [3, 1], isf_square_dict)  # 可测真实beta

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1])],
                           "cols": [[3, 2], [2, 2, 1]],
                           "isf": sp.Matrix([[-Ra(5) / 8, Ra(3) / 8], [-Ra(3) / 8, -Ra(5) / 8]])}
        self.isf_22 = (5, [3, 2], [3, 1, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1]), ([2, 2], [3, 1]), ([2, 2], [2, 1, 1])],
                           "cols": [([3, 1, 1], 1), ([3, 1, 1], 2), [2, 2, 1], [2, 1, 1, 1]],
                           "isf": sp.Matrix([[0, Ra(1) / 2, -Ra(1) / 3, Ra(1) / 6],
                                             [-Ra(3) / 5, Ra(1) / 10, 0, -Ra(3) / 10],
                                             [Ra(1) / 4, 0, -Ra(1) / 4, -Ra(1) / 2],
                                             [Ra(3) / 20, Ra(2) / 5, Ra(5) / 12, -Ra(1) / 30]])}
        self.isf_23 = (5, [3, 2], [3, 1, 1], [2, 1, 1], isf_square_dict)  # 可测真实beta  # 还能测first_no_0

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1])],
                           "cols": [[2, 1, 1, 1]],
                           "isf": sp.Matrix([[1]])}
        self.isf_24 = (5, [3, 2], [3, 1, 1], [1, 1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1]), ([2, 2], [2, 2])],
                           "cols": [[2, 1, 1, 1], [1, 1, 1, 1, 1]],
                           "isf": sp.Matrix([[-Ra(2) / 5, Ra(3) / 5], [Ra(3) / 5, Ra(2) / 5]])}
        self.isf_25 = (5, [3, 2], [2, 2, 1], [1, 1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1]), ([2, 2], [1, 1, 1, 1])],
                           "cols": [[3, 2], [2, 2, 1]],
                           "isf": sp.Matrix([[Ra(3) / 8, Ra(5) / 8], [Ra(5) / 8, -Ra(3) / 8]])}
        self.isf_26 = (5, [3, 2], [2, 1, 1, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [[4, 1], ([3, 2], 1), ([3, 2], 2), [3, 1, 1]],
                           "isf": sp.Matrix([[Ra(5) / 12, Ra(1) / 2, Ra(1) / 12, 0],
                                             [-Ra(1) / 12, 0, Ra(5) / 12, Ra(1) / 2],
                                             [-Ra(1) / 12, 0, Ra(5) / 12, -Ra(1) / 2],
                                             [Ra(5) / 12, -Ra(1) / 2, Ra(1) / 12, 0]])}
        self.isf_27 = (5, [3, 1, 1], [3, 1, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [([3, 2], 1), ([3, 2], 2), ([2, 2, 1], 1), ([2, 2, 1], 2)],
                           "isf": sp.Matrix([[-Ra(3) / 16, Ra(1) / 2, Ra(5) / 16, 0],
                                             [Ra(5) / 16, 0, Ra(3) / 16, Ra(1) / 2],
                                             [Ra(5) / 16, 0, Ra(3) / 16, -Ra(1) / 2],
                                             [Ra(3) / 16, Ra(1) / 2, -Ra(5) / 16, 0]])}
        self.isf_28 = (5, [3, 1, 1], [3, 1, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 1, 1]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [[3, 1, 1], ([2, 2, 1], 1), ([2, 2, 1], 2), [2, 1, 1, 1]],
                           "isf": sp.Matrix([[Ra(1) / 2, 0, -Ra(5) / 12, Ra(1) / 12],
                                             [0, -Ra(1) / 2, Ra(1) / 12, Ra(5) / 12],
                                             [0, -Ra(1) / 2, -Ra(1) / 12, -Ra(5) / 12],
                                             [Ra(1) / 2, 0, Ra(5) / 12, -Ra(1) / 12]])}
        self.isf_29 = (5, [3, 1, 1], [3, 1, 1], [2, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 2]), ([3, 1], [2, 1, 1]), ([2, 1, 1], [2, 2]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [[4, 1], [3, 2], ([3, 1, 1], 1), ([3, 1, 1], 2)],
                           "isf": sp.Matrix([[Ra(1) / 2, Ra(1) / 4, Ra(1) / 4, 0],
                                             [Ra(1) / 6, -Ra(1) / 3, 0, Ra(1) / 2],
                                             [Ra(1) / 30, -Ra(5) / 12, Ra(3) / 20, -Ra(2) / 5],
                                             [Ra(3) / 10, 0, -Ra(3) / 5, -Ra(1) / 10]])}
        self.isf_30 = (5, [3, 1, 1], [2, 2, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [[3, 2], [2, 2, 1]],
                           "isf": sp.Matrix([[-Ra(3) / 8, Ra(5) / 8], [Ra(5) / 8, Ra(3) / 8]])}
        self.isf_31 = (5, [3, 1, 1], [2, 2, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 2]), ([3, 1], [2, 1, 1]), ([2, 1, 1], [2, 2]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [([3, 1, 1], 1), ([3, 1, 1], 2), [2, 2, 1], [2, 1, 1, 1]],
                           "isf": sp.Matrix([[Ra(3) / 20, Ra(2) / 5, -Ra(5) / 12, Ra(1) / 30],
                                             [-Ra(3) / 5, Ra(1) / 10, 0, Ra(3) / 10],
                                             [-Ra(1) / 4, 0, -Ra(1) / 4, -Ra(1) / 2],
                                             [0, -Ra(1) / 2, -Ra(1) / 3, Ra(1) / 6]])}
        self.isf_32 = (5, [3, 1, 1], [2, 2, 1], [2, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [2, 1, 1]), ([2, 1, 1], [2, 1, 1]), ([2, 1, 1], [1, 1, 1, 1])],
                           "cols": [[4, 1], [3, 2], [3, 1, 1]],
                           "isf": sp.Matrix([[-Ra(2) / 3, Ra(5) / 24, Ra(1) / 8],
                                             [0, -Ra(3) / 8, Ra(5) / 8],
                                             [Ra(1) / 3, Ra(5) / 12, Ra(1) / 4]])}
        self.isf_33 = (5, [3, 1, 1], [2, 1, 1, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 2], [2, 1, 1]), ([2, 1, 1], [2, 2]), ([2, 1, 1], [2, 1, 1])],
                           "cols": [[4, 1], [3, 2], [3, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 3, Ra(1) / 6, Ra(1) / 2],
                                             [-Ra(1) / 3, Ra(1) / 6, -Ra(1) / 2],
                                             [Ra(1) / 3, Ra(2) / 3, 0]])}
        self.isf_34 = (5, [2, 2, 1], [2, 2, 1], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 2], [2, 2])], "cols": [[2, 1, 1, 1]], "isf": sp.Matrix([[-1]])}
        self.isf_35 = (5, [2, 2, 1], [2, 2, 1], [1, 1, 1, 1], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1, 1], [2, 1, 1])], "cols": [[4, 1]], "isf": sp.Matrix([[1]])}
        self.isf_36 = (5, [2, 2, 1], [2, 1, 1, 1], [4], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1, 1], [2, 1, 1]), ([1, 1, 1, 1], [1, 1, 1, 1])],
                           "cols": [[5], [4, 1]],
                           "isf": sp.Matrix([[Ra(3) / 4, Ra(1) / 4], [Ra(1) / 4, -Ra(3) / 4]])}
        self.isf_37 = (5, [2, 1, 1, 1], [2, 1, 1, 1], [4], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1, 1], [2, 1, 1])], "cols": [[3, 2]], "isf": sp.Matrix([[1]])}
        self.isf_38 = (5, [2, 1, 1, 1], [2, 1, 1, 1], [2, 2], isf_square_dict)

        isf_square_dict = {"rows": [([2, 1, 1], [2, 1, 1])], "cols": [[3, 1, 1]], "isf": sp.Matrix([[1]])}
        self.isf_39 = (5, [2, 1, 1, 1], [2, 1, 1, 1], [2, 1, 1], isf_square_dict)

        # 自己算的
        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 2]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 2])],
                           "cols": [[4, 1], [3, 2], ([3, 1, 1], 1), ([3, 1, 1], 2)],
                           "isf": sp.Matrix([[-Ra(3) / 10, 0, Ra(3) / 5, Ra(1) / 10],
                                             [-Ra(1) / 30, Ra(5) / 12, -Ra(3) / 20, Ra(2) / 5],
                                             [-Ra(1) / 6, Ra(1) / 3, 0, -Ra(1) / 2],
                                             [Ra(1) / 2, Ra(1) / 4, Ra(1) / 4, 0]])}
        self.isf_40 = (5, [3, 1, 1], [3, 2], [3, 1], isf_square_dict)

        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 2]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 2])],
                           "cols": [([3, 1, 1], 1), ([3, 1, 1], 2), [2, 2, 1], [2, 1, 1, 1]],
                           "isf": sp.Matrix([[0, Ra(1) / 2, -Ra(1) / 3, Ra(1) / 6],
                                             [Ra(1) / 4, 0, -Ra(1) / 4, -Ra(1) / 2],
                                             [-Ra(3) / 5, Ra(1) / 10, 0, -Ra(3) / 10],
                                             [Ra(3) / 20, Ra(2) / 5, Ra(5) / 12, -Ra(1) / 30]])}
        self.isf_41 = (5, [3, 1, 1], [3, 2], [2, 1, 1], isf_square_dict)

        # TODO 错误原因待定，目前看是因为isf_matrix顺序变了导致高斯消元法在计算有自由度的基失时，会有'旋转'，不再对称
        # 但，它依然是自正交归一的！
        # 它不是错了，只是不符合我们约定的对称性
        isf_square_dict = {"rows": [([3, 1], [3, 1]), ([3, 1], [2, 2]), ([2, 1, 1], [3, 1]), ([2, 1, 1], [2, 2])],
                           "cols": [([3, 1, 1], 1), ([3, 1, 1], 2), [2, 2, 1], [2, 1, 1, 1]],
                           "isf": sp.Matrix([[-Ra(1) / 14, -Ra(3) / 7, Ra(1) / 3, Ra(1) / 6],
                                             [-Ra(3) / 14, Ra(1) / 28, -Ra(1) / 4, Ra(1) / 2],
                                             [Ra(7) / 10, 0, 0, Ra(3) / 10],
                                             [-Ra(1) / 70, Ra(15) / 28, Ra(5) / 12, Ra(1) / 30]])}
        self.isf_41_wrong = (5, [3, 1, 1], [3, 2], [2, 1, 1], isf_square_dict)

        # cgc data
        # 格式 Sn, σ, μ, ν, β, m, cgc_square_dict
        cgc_square_dict = {(1, 1,): 1, "N": 1}
        self.cgc_1 = (1, [1], [1], [1], None, 1, cgc_square_dict)

        # 数据来源见《群表示论的新途径》陈金全（上海科学技术出版社1984）表4.13
        cgc_square_dict = {(1, 1,): 1, (2, 2): 1, "N": 2}
        self.cgc_2 = (3, [2, 1], [2, 1], [3], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1,): 1, (2, 2): -1, "N": 2}
        self.cgc_3 = (3, [2, 1], [2, 1], [2, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2,): -1, (2, 1): -1, "N": 2}
        self.cgc_4 = (3, [2, 1], [2, 1], [2, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2,): 1, (2, 1): -1, "N": 2}
        self.cgc_5 = (3, [2, 1], [2, 1], [1, 1, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1,): 1, (2, 2): 1, (3, 3): 1, "N": 3}
        self.cgc_6 = (4, [3, 1], [3, 1], [4], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1,): 4, (2, 2): -1, (3, 3): -1, "N": 6}
        self.cgc_7 = (4, [3, 1], [3, 1], [3, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2,): -1, (2, 1): -1, (2, 2): 2, (3, 3): -2, "N": 6}
        self.cgc_8 = (4, [3, 1], [3, 1], [3, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 3,): -1, (3, 1): -1, (2, 3): -2, (3, 2): -2, "N": 6}
        self.cgc_9 = (4, [3, 1], [3, 1], [3, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(2, 3,): 1, (3, 2): -1, "N": 2}
        self.cgc_10 = (4, [3, 1], [3, 1], [2, 1, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(2, 1): 1, (3, 2): 1, "N": 2}
        self.cgc_11 = (4, [3, 1], [2, 2], [3, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 2, (2, 1): 1, (3, 2): -1, "N": 4}
        self.cgc_12 = (4, [3, 1], [2, 2], [3, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): 2, (2, 2): -1, (3, 1): -1, "N": 4}
        self.cgc_13 = (4, [3, 1], [2, 2], [3, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(2, 2): -1, (3, 1): 1, "N": 2}
        self.cgc_14 = (4, [3, 1], [2, 2], [2, 1, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 1): 1, (2, 2): 1, "N": 2}
        self.cgc_15 = (4, [2, 2], [2, 2], [4], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 1, (2, 2): -1, "N": 2}
        self.cgc_16 = (4, [2, 2], [2, 2], [2, 2], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -1, (2, 1): -1, "N": 2}
        self.cgc_17 = (4, [2, 2], [2, 2], [2, 2], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): 1, (2, 1): -1, "N": 2}
        self.cgc_18 = (4, [2, 2], [2, 2], [1, 1, 1, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 1, (2, 2): 1, (3, 3): 1, (4, 4): 1, "N": 4}
        self.cgc_19 = (5, [4, 1], [4, 1], [5], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -3, (2, 1): -3, (2, 2): 20, (3, 3): -5, (4, 4): -5, "N": 36}
        self.cgc_20 = (5, [4, 1], [4, 1], [4, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): 15, (2, 1): 15, (2, 2): 4, (3, 3): -1, (4, 4): -1, "N": 36}
        self.cgc_21 = (5, [4, 1], [4, 1], [3, 2], None, 1, cgc_square_dict)

        cgc_square_dict = {(2, 3): 1, (3, 2): -1, "N": 2}
        self.cgc_22 = (5, [4, 1], [4, 1], [3, 1, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 4): -3, (2, 4): -5, (3, 4): -10, (4, 1): -3, (4, 2): -5, (4, 3): -10, "N": 36}
        self.cgc_23 = (5, [4, 1], [4, 1], [4, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(2, 4): 2, (3, 4): -1, (4, 2): 2, (4, 3): -1, "N": 6}
        self.cgc_24 = (5, [4, 1], [4, 1], [3, 2], None, 5, cgc_square_dict)

        cgc_square_dict = {(2, 4): 1, (4, 2): -1, "N": 2}
        self.cgc_25 = (5, [4, 1], [4, 1], [3, 1, 1], None, 5, cgc_square_dict)

        cgc_square_dict = {(1, 2): 15, (2, 2): -1,
                           (2, 4): 12, (3, 1): -1, (3, 2): 2, (3, 4): 6, (4, 3): -2, (4, 5): -6, "N": 45}
        self.cgc_26 = (5, [4, 1], [3, 2], [4, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 1): 12, (2, 1): 20, (3, 2): -5, (3, 4): -15, (4, 3): -5, (4, 5): -15, "N": 72}
        self.cgc_27 = (5, [4, 1], [3, 2], [3, 2], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 20, (2, 1): -12, (3, 2): 3, (3, 4): -1, (4, 3): 3, (4, 5): -1, "N": 40}
        self.cgc_28 = (5, [4, 1], [3, 2], [3, 1, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(2, 2): 6, (2, 4): 2, (3, 1): -6, (3, 4): -1, (4, 5): 1, "N": 16}
        self.cgc_29 = (5, [4, 1], [3, 2], [2, 2, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): 15, (2, 3): -1, (2, 5): 12, (3, 3): -2,
                           (3, 5): -6, (4, 1): -1, (4, 2): -2, (4, 4): -6, "N": 45}
        self.cgc_30 = (5, [4, 1], [3, 2], [4, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 3): 40, (2, 3): 6, (2, 5): -2, (3, 3): 12,
                           (3, 5): 1, (4, 1): 6, (4, 2): 12, (4, 4): 1, "N": 80}
        self.cgc_31 = (5, [4, 1], [3, 2], [3, 1, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(3, 3): 3, (3, 5): -1, (4, 2): -3, (4, 4): 1, "N": 8}
        self.cgc_32 = (5, [4, 1], [3, 2], [2, 2, 1], None, 5, cgc_square_dict)

        cgc_square_dict = {(1, 1): -1, (4, 3): 1, (5, 4): 1, "N": 3}
        self.cgc_33 = (5, [3, 1, 1], [4, 1], [4, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(1, 1): 20, (1, 2): -12, (2, 3): 3, (3, 4): 3, (4, 3): 5, (5, 4): 5, "N": 48}
        self.cgc_34 = (5, [3, 1, 1], [4, 1], [3, 2], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 3): 3, (2, 2): -3, (4, 1): -12, (4, 2): 5, (4, 3): 10, (5, 4): -10, (6, 4): 5, "N": 48}
        self.cgc_35 = (5, [3, 1, 1], [4, 1], [3, 1, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 3): 30, (2, 2): 30,
                           (2, 3): 15, (3, 4): -15, (4, 2): 2, (4, 3): -1, (5, 4): 1, (6, 4): 2, "N": 96}
        self.cgc_36 = (5, [3, 1, 1], [4, 1], [2, 2, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(3, 1): -1, (5, 2): -1, (6, 3): -1, "N": 3}
        self.cgc_37 = (5, [3, 1, 1], [4, 1], [4, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 4): 3, (3, 2): -3,
                           (4, 4): -10, (5, 1): -12, (5, 2): 5, (5, 3): -10, (6, 3): -5, "N": 48}
        self.cgc_38 = (5, [3, 1, 1], [4, 1], [3, 1, 1], None, 5, cgc_square_dict)

        cgc_square_dict = {(1, 4): 30, (2, 4): -15, (3, 2): 30, (3, 3): -15,
                           (4, 4): 1, (5, 2): 2, (5, 3): 1, (6, 3): -2, "N": 96}
        self.cgc_39 = (5, [3, 1, 1], [4, 1], [2, 2, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(4, 4): 1, (5, 3): -1, (6, 2): 1, "N": 3}
        self.cgc_40 = (5, [3, 1, 1], [4, 1], [2, 1, 1, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 2): -2, (1, 4): -6, (2, 1): -2, (2, 2): 4,
                           (2, 4): -3, (3, 3): -4, (3, 5): 3, (4, 1): -6, (4, 2): -3, (5, 3): 3, "N": 36}
        self.cgc_41 = (5, [3, 2], [3, 2], [4, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 2): -8, (1, 4): 6, (2, 1): -8, (2, 2): 16,
                           (2, 4): 3, (3, 3): -16, (3, 5): -3, (4, 1): 6, (4, 2): 3, (5, 3): -3, "N": 72}
        self.cgc_42 = (5, [3, 2], [3, 2], [3, 2], None, 2, cgc_square_dict)

        cgc_square_dict = {(2, 4): 1, (3, 5): 1, (4, 2): -1, (5, 3): -1, "N": 4}
        self.cgc_43 = (5, [3, 2], [3, 2], [3, 1, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 4): 2, (2, 4): -1, (3, 5): 1, (4, 1): 2, (4, 2): -1, (5, 3): 1, "N": 8}
        self.cgc_44 = (5, [3, 2], [3, 2], [2, 2, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): -2, (1, 5): -6, (2, 3): -4, (2, 5): 3, (3, 1): -2,
                           (3, 2): -4, (3, 4): 3, (4, 3): 3, (5, 1): -6, (5, 2): 3, "N": 36}
        self.cgc_45 = (5, [3, 2], [3, 2], [4, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 3): 8, (1, 5): 6, (2, 5): 3, (3, 1): -8,
                           (3, 4): 3, (4, 3): -3, (5, 1): -6, (5, 2): -3, "N": 40}
        self.cgc_46 = (5, [3, 2], [3, 2], [3, 1, 1], None, 5, cgc_square_dict)

        cgc_square_dict = {(2, 5): -1, (3, 4): 1, (4, 3): 1, (5, 2): -1, "N": 4}
        self.cgc_47 = (5, [3, 2], [3, 2], [2, 2, 1], None, 5, cgc_square_dict)

        cgc_square_dict = {(4, 5): 1, (5, 4): -1, "N": 2}
        self.cgc_48 = (5, [3, 2], [3, 2], [2, 1, 1, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 2): 6, (1, 4): -2, (2, 1): 6, (2, 2): -12, (2, 4): -1, (3, 3): 12, (3, 5): 1, (4, 1): 10,
                           (4, 4): -15, (5, 5): 15, (6, 3): -10, (6, 5): -30, "N": 120}
        self.cgc_49 = (5, [3, 1, 1], [3, 2], [4, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(2, 4): 5, (3, 5): 5, (4, 2): 4, (4, 4): 3, (5, 3): 4, (5, 5): 3, "N": 24}
        self.cgc_50 = (5, [3, 1, 1], [3, 2], [3, 2], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 16, (2, 2): -4, (2, 4): -3, (3, 3): -4, (3, 5): -3, (4, 4): 5, (5, 5): 5, "N": 40}
        self.cgc_51 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 1, cgc_square_dict)  # 可以测beta

        cgc_square_dict = {(1, 2): -8, (1, 4): -6, (2, 1): -8, (2, 2): 16, (2, 4): -3, (3, 3): -16, (3, 5): 3,
                           (4, 4): -5, (5, 5): 5, (6, 5): -10, "N": 80}
        self.cgc_52 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 2, cgc_square_dict)

        cgc_square_dict = {(1, 4): 10, (2, 4): -5, (3, 5): 5, (4, 1): -8, (4, 2): -16,
                           (4, 4): 3, (5, 3): 16, (5, 5): -3, (6, 3): -8, (6, 5): -6, "N": 80}
        self.cgc_53 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 4, cgc_square_dict)

        cgc_square_dict = {(1, 1): 4, (2, 2): -1, (2, 4): 12, (3, 3): -1, (3, 5): 12, (4, 2): -15, (5, 3): -15, "N": 60}
        self.cgc_54 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -1, (1, 4): 12, (2, 1): -1, (2, 2): 2,
                           (2, 4): 6, (3, 3): -2, (3, 5): -6, (4, 1): 15, (6, 3): -15, "N": 60}
        self.cgc_55 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): -15, (2, 1): 15, (4, 1): 1, (4, 2): 2,
                           (4, 4): 6, (5, 3): -2, (5, 5): -6, (6, 3): 1, (6, 5): -12, "N": 60}
        self.cgc_56 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 4, cgc_square_dict)

        cgc_square_dict = {(1, 2): 8, (1, 4): -6, (2, 1): -8, (2, 4): 3, (3, 5): -3,
                           (4, 4): 5, (5, 5): -5, (6, 5): -10, "N": 48}
        self.cgc_57 = (5, [3, 1, 1], [3, 2], [2, 2, 1], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 2): 10, (1, 4): 30, (2, 1): -10, (2, 4): -15, (3, 5): 15, (4, 1): 6, (4, 2): 12,
                           (4, 4): 1, (5, 3): -12, (5, 5): -1, (6, 3): 6, (6, 5): -2, "N": 120}
        self.cgc_58 = (5, [3, 1, 1], [3, 2], [2, 1, 1, 1], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 5): 10, (2, 5): -5, (3, 4): -5, (4, 5): 3, (5, 1): -8,
                           (5, 4): 3, (6, 2): -8, (6, 4): 6, "N": 48}
        self.cgc_59 = (5, [3, 1, 1], [3, 2], [3, 2], None, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): -8, (1, 5): -6, (2, 3): -16, (2, 5): 3, (3, 1): -8,
                           (3, 2): -16, (3, 4): 3, (4, 5): 5, (5, 4): 5, (6, 4): 10, "N": 80}
        self.cgc_60 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 3, cgc_square_dict)

        cgc_square_dict = {(1, 5): 10, (2, 5): 5, (3, 4): 5, (4, 3): 16, (4, 5): -3, (5, 1): -8,
                           (5, 2): 16, (5, 4): -3, (6, 2): 8, (6, 4): 6, "N": 80}
        self.cgc_61 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 5, cgc_square_dict)

        cgc_square_dict = {(2, 5): -5, (3, 4): 5, (4, 3): -4, (4, 5): -3, (5, 2): 4, (5, 4): 3, (6, 1): 16, "N": 40}
        self.cgc_62 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 1, 6, cgc_square_dict)

        cgc_square_dict = {(1, 3): -1, (1, 5): 12, (2, 3): -2, (2, 5): -6, (3, 1): -1,
                           (3, 2): -2, (3, 4): -6, (5, 1): 15, (6, 2): 15, "N": 60}
        self.cgc_63 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): -15, (3, 1): 15, (4, 3): -2, (4, 5): -6, (5, 1): 1,
                           (5, 2): -2, (5, 4): -6, (5, 1): 1, (5, 2): -2, (5, 4): -6, (6, 2): -1, (6, 4): 12, "N": 60}
        self.cgc_64 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 5, cgc_square_dict)

        cgc_square_dict = {(2, 3): -15, (3, 2): 15, (4, 3): 1, (4, 5): -12,
                           (5, 2): -1, (5, 4): 12, (6, 1): -4, "N": 60}
        self.cgc_65 = (5, [3, 1, 1], [3, 2], [3, 1, 1], 2, 6, cgc_square_dict)

        cgc_square_dict = {(1, 3): 6, (2, 3): -3, (3, 1): 6,
                           (3, 2): -3, (4, 3): -5, (5, 1): -10, (5, 2): -5, (6, 2): 10, "N": 48}
        self.cgc_66 = (5, [3, 1, 1], [3, 2], [2, 2, 1], None, 2, cgc_square_dict)

        cgc_square_dict = {(4, 3): -1, (5, 2): 1, (6, 1): -1, "N": 3}
        self.cgc_67 = (5, [3, 1, 1], [3, 2], [2, 1, 1, 1], None, 4, cgc_square_dict)

        cgc_square_dict = {(1, 1): 1, (2, 2): 1, (3, 3): 1, (4, 4): 1, (5, 5): 1, (6, 6): 1, "N": 6}
        self.cgc_68 = (5, [3, 1, 1], [3, 1, 1], [5], None, 1, cgc_square_dict)

        cgc_square_dict = {(1, 1): 4, (2, 2): -1, (3, 3): -1, (4, 4): -1, (5, 5): -1, (6, 6): 4, "N": 12}
        self.cgc_69 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 1, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -1, (2, 1): -1, (2, 2): 2, (3, 3): -2,
                           (4, 4): -2, (5, 5): 2, (5, 6): 1, (6, 5): 1, "N": 12}
        self.cgc_70 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 1, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): -6, (1, 4): 10, (2, 1): -6, (2, 2): -3,
                           (2, 4): -5, (3, 3): 3, (3, 5): 5, (3, 6): 10, (4, 1): 10,
                           (4, 2): -5, (4, 4): 3, (5, 3): 5, (5, 5): -3, (5, 6): 6, (6, 3): 10, (6, 5): 6, "N": 96}
        self.cgc_71 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 1, 4, cgc_square_dict)

        cgc_square_dict = {(1, 1): 4, (2, 2): -1, (2, 4): 15, (3, 3): -1, (3, 5): 15,
                           (4, 2): 15, (4, 4): 1, (5, 3): 15, (5, 5): 1, (6, 6): -4, "N": 72}
        self.cgc_72 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 2, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -1, (1, 4): -15, (2, 1): -1, (2, 2): 2, (3, 3): -2, (3, 6): 15, (4, 1): -15,
                           (4, 4): 2, (5, 5): -2, (5, 6): -1, (6, 3): 15, (6, 5): -1, "N": 72}
        self.cgc_73 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 2, 2, cgc_square_dict)

        cgc_square_dict = {(1, 2): 2, (2, 1): 2, (2, 2): 1, (3, 3): -1,
                           (4, 4): 1, (5, 5): -1, (5, 6): 2, (6, 5): 2, "N": 12}
        self.cgc_74 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 2, 4, cgc_square_dict)

        cgc_square_dict = {(1, 2): 10, (1, 4): 6, (2, 1): 10, (2, 2): 5,
                           (2, 4): -3, (3, 3): -5, (3, 5): 3, (3, 6): 6, (4, 1): 6,
                           (4, 2): -3, (4, 4): -5, (5, 3): 3, (5, 5): 5, (5, 6): -10,
                           (6, 3): 6, (6, 5): -10, "N": 96}
        self.cgc_75 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 1, 1, cgc_square_dict)

        cgc_square_dict = {(1, 4): -1, (2, 4): -2, (3, 5): 2, (3, 6): -1, (4, 1): -1,
                           (4, 2): -2, (5, 3): 2, (6, 3): -1, "N": 12}
        self.cgc_76 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 1, 3, cgc_square_dict)

        cgc_square_dict = {(1, 4): 2, (2, 4): -1, (3, 5): 1, (3, 6): 2, (4, 1): -2,
                           (4, 2): 1, (5, 3): -1, (6, 3): -2, "N": 12}
        self.cgc_77 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 2, 1, cgc_square_dict)

        cgc_square_dict = {(1, 2): -15, (1, 4): 1, (2, 1): 15, (2, 4): 2, (3, 5): -2, (3, 6): 1, (4, 1): -1,
                           (4, 2): -2, (5, 3): 2, (5, 6): 15, (6, 3): -1, (6, 5): -15, "N": 72}
        self.cgc_78 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 2, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): -1, (2, 3): -2, (3, 1): -1, (3, 2): -2,
                           (4, 5): 2, (4, 6): -1, (5, 4): 2, (6, 4): -1, "N": 12}
        self.cgc_79 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 1, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): -6, (1, 5): 10, (2, 3): 3, (2, 5): 5,
                           (2, 6): -10, (3, 1): -6, (3, 2): 3, (3, 4): 5, (4, 3): 5,
                           (4, 5): -3, (4, 6): -6, (5, 1): 10, (5, 2): 5, (5, 4): -3,
                           (6, 2): -10, (6, 4): -6, "N": 96}
        self.cgc_80 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 1, 5, cgc_square_dict)

        cgc_square_dict = {(1, 3): -1, (1, 5): -15, (2, 3): -2, (2, 6): -15, (3, 1): -1, (3, 2): -2,
                           (4, 5): -2, (4, 6): 1, (5, 1): -15, (5, 4): -2, (6, 2): -15, (6, 4): 1, "N": 72}
        self.cgc_81 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 2, 3, cgc_square_dict)

        cgc_square_dict = {(1, 3): 2, (2, 3): -1, (3, 1): 2, (3, 2): -1,
                           (4, 5): -1, (4, 6): -2, (5, 4): -1, (6, 4): -2, "N": 12}
        self.cgc_82 = (5, [3, 1, 1], [3, 1, 1], [3, 2], 2, 5, cgc_square_dict)

        cgc_square_dict = {(1, 3): 10, (1, 5): 6, (2, 3): -5, (2, 5): 3,
                           (2, 6): -6, (3, 1): 10, (3, 2): -5, (3, 4): 3, (4, 3): 3,
                           (4, 5): 5, (4, 6): 10, (5, 1): 6, (5, 2): 3, (5, 4): 5, (6, 2): -6, (6, 4): 10, "N": 96}
        self.cgc_83 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 1, 2, cgc_square_dict)

        cgc_square_dict = {(1, 5): -1, (2, 5): 2, (2, 6): 1, (3, 4): 2, (4, 3): 2,
                           (5, 1): -1, (5, 2): 2, (6, 2): 1, "N": 12}
        self.cgc_84 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 1, 4, cgc_square_dict)

        cgc_square_dict = {(1, 6): 4, (2, 5): 1, (3, 4): -1, (4, 3): -1, (5, 2): 1, (6, 1): 4, "N": 12}
        self.cgc_85 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 1, 5, cgc_square_dict)

        cgc_square_dict = {(1, 5): 2, (2, 5): 1, (2, 6): -2, (3, 4): 1, (4, 3): -1,
                           (5, 1): -2, (5, 2): -1, (6, 2): 2, "N": 12}
        self.cgc_86 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 2, 2, cgc_square_dict)

        cgc_square_dict = {(1, 3): -15, (1, 5): 1, (2, 5): -2, (2, 6): -1, (3, 1): 15, (3, 4): -2, (4, 3): 2,
                           (4, 6): -15, (5, 1): -1, (5, 2): 2, (6, 2): 1, (6, 4): 15, "N": 72}
        self.cgc_87 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 2, 4, cgc_square_dict)

        cgc_square_dict = {(1, 6): -4, (2, 3): -15, (2, 5): -1, (3, 2): 15, (3, 4): 1, (4, 3): -1,
                           (4, 5): 15, (5, 2): 1, (5, 4): -15, (6, 1): 4, "N": 72}
        self.cgc_88 = (5, [3, 1, 1], [3, 1, 1], [2, 2, 1], 2, 5, cgc_square_dict)

        # isf matrix
        # 格式 Sn, σ, μ, ν_st, row_index_tmp_list, isf_matrix
        row_index_tmp_list = [([2], [2], None), ([1, 1], [1, 1], None)]
        isf_matrix = sp.Matrix([[Ra(1) / 2, Ra(3) / 2], [Ra(3) / 2, Ra(1) / 2]])
        self.isf_matrix_1 = (3, [2, 1], [2, 1], [2], row_index_tmp_list, isf_matrix)

        row_index_tmp_list = [([3], [2, 1], None), ([2, 1], [2, 1], None)]
        isf_matrix = sp.Matrix([[0, 2], [2, 0]])
        self.isf_matrix_2 = (4, [3, 1], [2, 2], [2, 1], row_index_tmp_list, isf_matrix)  # 4-19 例2


class TestISFAndCGC(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_isf_and_cgc_protected")
        self.protector.protector_setup()

        # 准备前文
        s_n = 5
        self.test_sn = s_n
        flag, msg = create_young_diagrams(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_branching_laws(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_young_tableaux(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_yamanouchi_matrix(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_characters_and_gi(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_cg_series(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_eigenvalues(s_n)
        assert flag
        assert msg == s_n

        self.decimals = len(str(isf_0_error_value))

        self.data = Data()
        self.helper = Helper()

        # self.isf_ban_set = {21, 23, 27, 28, 29, 30, 32}
        # self.isf_ban_set = {30, 32}
        self.isf_ban_set = set()
        self.isf_num_list = list(set(range(1, 39 + 1)) - self.isf_ban_set)

        # self.cgc_ban_set = set(chain(range(51, 56 + 1), range(60, 65 + 1), range(69, 88 + 1)))
        self.cgc_ban_set = set()
        self.cgc_num_list = list(set(range(1, 88 + 1)) - self.cgc_ban_set)
        # self.cgc_num_list = list(range(1, 88 + 1))

        _, self.isf_s_n_finish_file_name = get_isf_finish_s_n_name()
        _, self.isf_s_n_finish_full_file_name = get_isf_finish_s_n_name(is_full_path=True)
        self.isf_create_time_dict = {}  # 用于检查计算好的部分不会重复计算

        _, self.cgc_s_n_finish_file_name = get_cgc_finish_s_n_name()
        _, self.cgc_s_n_finish_full_file_name = get_cgc_finish_s_n_name(is_full_path=True)
        self.cgc_create_time_dict = {}  # 用于检查计算好的部分不会重复计算

        self.isf_matrix_num_list = list(range(1, 2 + 1))

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    # @pytest.mark.skip("pass")
    def test_001_create_isf_and_cgc_s_n_1(self):  # 初始cgc
        """for s_n=1, there is no finish db"""
        # check with no db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.isf_s_n_finish_full_file_name + ex)
            assert not os.path.exists(self.cgc_s_n_finish_full_file_name + ex)
        flag, isf_square_dict = load_isf(*self.data.isf_1[: -1], is_flag_true_if_not_s_n=True)
        assert flag
        assert isf_square_dict is False
        flag, isf_square_dict = load_isf(*self.data.isf_1[: -1], is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(isf_square_dict, str)
        flag, cgc_square_dict = load_cgc(*self.data.cgc_1[: -1], is_flag_true_if_not_s_n=True)
        assert flag
        assert cgc_square_dict is False
        flag, cgc_square_dict = load_cgc(*self.data.cgc_1[: -1], is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(cgc_square_dict, str)

        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == 1
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == 0

        # check create_isf_and_cgc_s_n_1
        flag, finish_s_n = create_isf_and_cgc(1)
        assert flag
        assert finish_s_n == 1

        # check answer
        for nb in self.isf_num_list:
            isf_tuple = eval("self.data.isf_{}".format(nb))
            isf_param = isf_tuple[: -1]
            isf_answer = isf_tuple[-1]
            if isf_param[0] > finish_s_n:
                continue
            _, file_name = get_isf_file_name(*isf_param)
            flag, data = ISFInfo(isf_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_isf_file_name(*isf_param, is_full_path=True)
            _, full_finish_file_name = get_isf_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, isf = load_isf(*isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert isf["rows"] == isf_answer["rows"]
            assert isf["cols"] == isf_answer["cols"]
            assert isf["isf"] == isf_answer["isf"], "self.data.isf_{}, ={} with {}".format(nb, isf_param, isf_answer)

            # 以下都是和isf_x σ, μ对称isf的验证
            symmetry_isf_tuple = self.helper.calc_symmetry_isf(*isf_tuple)
            symmetry_isf_param = symmetry_isf_tuple[: -1]
            symmetry_isf_answer = symmetry_isf_tuple[-1]
            flag, symmetry_isf = load_isf(*symmetry_isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert symmetry_isf["rows"] == symmetry_isf_answer["rows"]
            assert symmetry_isf["cols"] == symmetry_isf_answer["cols"]
            assert symmetry_isf["isf"] == symmetry_isf_answer["isf"], \
                "self.data.symmetry_isf_{}, ={} with {}".format(nb, symmetry_isf_param, symmetry_isf_answer)

        for nb in self.cgc_num_list:
            cgc_tuple = eval("self.data.cgc_{}".format(nb))
            cgc_param = cgc_tuple[: -1]
            cgc_answer = cgc_tuple[-1]
            if cgc_param[0] > finish_s_n:
                continue
            _, file_name = get_cgc_file_name(*cgc_param)
            flag, data = CGCInfo(cgc_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            self.cgc_create_time_dict["S1"] = data.get("create_time")
            _, full_file_name = get_cgc_file_name(*cgc_param, is_full_path=True)
            _, full_finish_file_name = get_cgc_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cgc = load_cgc(*cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            cgc_answer_n = cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in cgc_answer.values()) == cgc_answer_n
            assert 1 - isf_0_error_value < abs(cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in cgc_answer.items():
                assert cgc[cgc_k] == Ra(cgc_v)/cgc_answer_n, \
                    "self.data.cgc_{}, ={} with {} with key={} and cgc={}".format(nb, cgc_param, cgc_answer, cgc_k, cgc)
            cgc_answer["N"] = cgc_answer_n

            # 以下都是和cgc_x σ, μ对称cgc的验证
            symmetry_cgc_tuple = self.helper.calc_symmetry_cgc(*cgc_tuple)
            symmetry_cgc_param = symmetry_cgc_tuple[: -1]
            symmetry_cgc_answer = symmetry_cgc_tuple[-1]
            flag, symmetry_cgc = load_cgc(*symmetry_cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            symmetry_cgc_answer_n = symmetry_cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in symmetry_cgc_answer.values()) == symmetry_cgc_answer_n
            assert 1 - isf_0_error_value < abs(symmetry_cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in symmetry_cgc_answer.items():
                assert symmetry_cgc[cgc_k] == Ra(cgc_v) / symmetry_cgc_answer_n, \
                    "symmetry_cgc_{}, ={} with {} with key={} and cgc={}".format(
                        nb, symmetry_cgc_param, symmetry_cgc_answer, cgc_k, cgc)
            symmetry_cgc_answer["N"] = symmetry_cgc_answer_n

        # check finish s_n
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == 1
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == 1

        _, cgc_finish_file_name = get_cgc_finish_s_n_name()
        flag, data = CGCInfo(cgc_finish_s_n).query_by_file_name(cgc_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        assert isinstance(data.get("flags").get("history_times").get("S1"), int)
        assert 0 <= data.get("flags").get("history_times").get("S1") <= 1
        flag, data_txt = CGCInfo(cgc_finish_s_n).query_txt_by_file_name(cgc_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        assert isinstance(data.get("history_times").get("S{}".format(1)), int)

    # @pytest.mark.skip("pass")
    def test_002_create_isf_and_cgc_s_n_2(self):  # 初始化isf
        # check with no isf db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.isf_s_n_finish_full_file_name + ex)
        flag, isf_square_dict = load_isf(*self.data.isf_1[: -1], is_flag_true_if_not_s_n=True)
        assert flag
        assert isf_square_dict is False
        flag, isf_square_dict = load_isf(*self.data.isf_1[: -1], is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(isf_square_dict, str)

        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == 1
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == 1

        # check create_isf_and_cgc_s_n_2
        flag, finish_s_n = create_isf_and_cgc(2)
        assert flag
        assert finish_s_n == 2

        # check create time
        _, file_name = get_cgc_file_name(*self.data.cgc_1[: -1])
        flag, data = CGCInfo(self.data.cgc_1[0]).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.cgc_create_time_dict["S1"] == data.get("create_time")

        # check answer
        for nb in self.isf_num_list:
            isf_tuple = eval("self.data.isf_{}".format(nb))
            isf_param = isf_tuple[: -1]
            isf_answer = isf_tuple[-1]
            if isf_param[0] > finish_s_n:
                continue
            _, file_name = get_isf_file_name(*isf_param)
            flag, data = ISFInfo(isf_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            self.isf_create_time_dict["S2"] = data.get("create_time")
            _, full_file_name = get_isf_file_name(*isf_param, is_full_path=True)
            _, full_finish_file_name = get_isf_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, isf = load_isf(*isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert isf["rows"] == isf_answer["rows"]
            assert isf["cols"] == isf_answer["cols"]
            assert isf["isf"] == isf_answer["isf"], "self.data.isf_{}, ={} with {}".format(nb, isf_param, isf_answer)

            # 以下都是和isf_x σ, μ对称isf的验证
            symmetry_isf_tuple = self.helper.calc_symmetry_isf(*isf_tuple)
            symmetry_isf_param = symmetry_isf_tuple[: -1]
            symmetry_isf_answer = symmetry_isf_tuple[-1]
            flag, symmetry_isf = load_isf(*symmetry_isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert symmetry_isf["rows"] == symmetry_isf_answer["rows"]
            assert symmetry_isf["cols"] == symmetry_isf_answer["cols"]
            assert symmetry_isf["isf"] == symmetry_isf_answer["isf"], \
                "self.data.symmetry_isf_{}, ={} with {}".format(nb, symmetry_isf_param, symmetry_isf_answer)

        for nb in self.cgc_num_list:
            cgc_tuple = eval("self.data.cgc_{}".format(nb))
            cgc_param = cgc_tuple[: -1]
            cgc_answer = cgc_tuple[-1]
            if cgc_param[0] > finish_s_n:
                continue
            _, file_name = get_cgc_file_name(*cgc_param)
            flag, data = CGCInfo(cgc_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_cgc_file_name(*cgc_param, is_full_path=True)
            _, full_finish_file_name = get_cgc_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cgc = load_cgc(*cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            cgc_answer_n = cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in cgc_answer.values()) == cgc_answer_n
            assert 1 - isf_0_error_value < abs(cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in cgc_answer.items():
                assert cgc[cgc_k] == Ra(cgc_v)/cgc_answer_n, \
                    "self.data.cgc_{}, ={} with {} with key={} and cgc={}".format(nb, cgc_param, cgc_answer, cgc_k, cgc)
            cgc_answer["N"] = cgc_answer_n

            # 以下都是和cgc_x σ, μ对称cgc的验证
            symmetry_cgc_tuple = self.helper.calc_symmetry_cgc(*cgc_tuple)
            symmetry_cgc_param = symmetry_cgc_tuple[: -1]
            symmetry_cgc_answer = symmetry_cgc_tuple[-1]
            flag, symmetry_cgc = load_cgc(*symmetry_cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            symmetry_cgc_answer_n = symmetry_cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in symmetry_cgc_answer.values()) == symmetry_cgc_answer_n
            assert 1 - isf_0_error_value < abs(symmetry_cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in symmetry_cgc_answer.items():
                assert symmetry_cgc[cgc_k] == Ra(cgc_v) / symmetry_cgc_answer_n, \
                    "symmetry_cgc_{}, ={} with {} with key={} and cgc={}".format(
                        nb, symmetry_cgc_param, symmetry_cgc_answer, cgc_k, cgc)
            symmetry_cgc_answer["N"] = symmetry_cgc_answer_n

        # check finish s_n
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == 2
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == 2

        # history_times
        _, isf_finish_file_name = get_isf_finish_s_n_name()
        flag, data = ISFInfo(isf_finish_s_n).query_by_file_name(isf_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = ISFInfo(isf_finish_s_n).query_txt_by_file_name(isf_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

        _, cgc_finish_file_name = get_cgc_finish_s_n_name()
        flag, data = CGCInfo(cgc_finish_s_n).query_by_file_name(cgc_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = CGCInfo(cgc_finish_s_n).query_txt_by_file_name(cgc_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    # @pytest.mark.skip("pass")
    def test_003_create_isf_and_cgc_s_n_3_to_4(self):
        start_sn = 3
        end_sn = 4
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == start_sn - 1
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == start_sn - 1

        # check create_isf_and_cgc_s_n_2
        flag, finish_s_n = create_isf_and_cgc(end_sn)
        assert flag
        assert finish_s_n == end_sn

        # check create time
        _, file_name = get_isf_file_name(*self.data.isf_1[: -1])
        flag, data = ISFInfo(self.data.isf_1[0]).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.isf_create_time_dict["S2"] == data.get("create_time")

        _, file_name = get_cgc_file_name(*self.data.cgc_1[: -1])
        flag, data = CGCInfo(self.data.cgc_1[0]).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.cgc_create_time_dict["S1"] == data.get("create_time")

        # check answer
        for nb in self.isf_num_list:
            isf_tuple = eval("self.data.isf_{}".format(nb))
            isf_param = isf_tuple[: -1]
            isf_answer = isf_tuple[-1]
            if isf_param[0] > end_sn or isf_param[0] < start_sn:
                continue
            _, file_name = get_isf_file_name(*isf_param)
            flag, data = ISFInfo(isf_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_isf_file_name(*isf_param, is_full_path=True)
            _, full_finish_file_name = get_isf_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, isf = load_isf(*isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert isf["rows"] == isf_answer["rows"]
            assert isf["cols"] == isf_answer["cols"]
            assert isf["isf"] == isf_answer["isf"], "self.data.isf_{}, ={} with {}".format(nb, isf_param, isf_answer)

            # 以下都是和isf_x σ, μ对称isf的验证
            symmetry_isf_tuple = self.helper.calc_symmetry_isf(*isf_tuple)
            symmetry_isf_param = symmetry_isf_tuple[: -1]
            symmetry_isf_answer = symmetry_isf_tuple[-1]
            flag, symmetry_isf = load_isf(*symmetry_isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert symmetry_isf["rows"] == symmetry_isf_answer["rows"]
            assert symmetry_isf["cols"] == symmetry_isf_answer["cols"]
            assert symmetry_isf["isf"] == symmetry_isf_answer["isf"], \
                "self.data.symmetry_isf_{}, ={} with {}".format(nb, symmetry_isf_param, symmetry_isf_answer)

        for nb in self.cgc_num_list:
            cgc_tuple = eval("self.data.cgc_{}".format(nb))
            cgc_param = cgc_tuple[: -1]
            cgc_answer = cgc_tuple[-1]
            if cgc_param[0] > end_sn or cgc_param[0] < start_sn:
                continue
            _, file_name = get_cgc_file_name(*cgc_param)
            flag, data = CGCInfo(cgc_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_cgc_file_name(*cgc_param, is_full_path=True)
            _, full_finish_file_name = get_cgc_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cgc = load_cgc(*cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            cgc_answer_n = cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in cgc_answer.values()) == cgc_answer_n
            assert cgc["N"] == 1
            for cgc_k, cgc_v in cgc_answer.items():
                assert cgc[cgc_k] == Ra(cgc_v) / cgc_answer_n, \
                    "self.data.cgc_{}, ={} with {} with key={} and cgc={}".format(nb, cgc_param, cgc_answer, cgc_k, cgc)
            cgc_answer["N"] = cgc_answer_n

            # 以下都是和cgc_x σ, μ对称cgc的验证
            symmetry_cgc_tuple = self.helper.calc_symmetry_cgc(*cgc_tuple)
            symmetry_cgc_param = symmetry_cgc_tuple[: -1]
            symmetry_cgc_answer = symmetry_cgc_tuple[-1]
            flag, symmetry_cgc = load_cgc(*symmetry_cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            symmetry_cgc_answer_n = symmetry_cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in symmetry_cgc_answer.values()) == symmetry_cgc_answer_n
            assert 1 - isf_0_error_value < abs(symmetry_cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in symmetry_cgc_answer.items():
                assert symmetry_cgc[cgc_k] == Ra(cgc_v) / symmetry_cgc_answer_n, \
                    "symmetry_cgc_{}, ={} with {} with key={} and cgc={}".format(
                        nb, symmetry_cgc_param, symmetry_cgc_answer, cgc_k, cgc)
            symmetry_cgc_answer["N"] = symmetry_cgc_answer_n

        # check finish s_n
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == end_sn
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == end_sn

        # history_times
        _, isf_finish_file_name = get_isf_finish_s_n_name()
        flag, data = ISFInfo(isf_finish_s_n).query_by_file_name(isf_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = ISFInfo(isf_finish_s_n).query_txt_by_file_name(isf_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

        _, cgc_finish_file_name = get_cgc_finish_s_n_name()
        flag, data = CGCInfo(cgc_finish_s_n).query_by_file_name(cgc_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        flag, data_txt = CGCInfo(cgc_finish_s_n).query_txt_by_file_name(cgc_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    # @pytest.mark.skip("pass")
    def test_004_create_isf_and_cgc_s_n_5_to_6(self):  # beta  # 005就该优化7～9了 放在benchmark里
        start_sn = 5
        end_sn = 5
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == start_sn - 1
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == start_sn - 1

        # check create_isf_and_cgc_s_n_2
        flag, finish_s_n = create_isf_and_cgc(end_sn)
        assert flag
        assert finish_s_n == end_sn

        # check create time
        _, file_name = get_isf_file_name(*self.data.isf_1[: -1])
        flag, data = ISFInfo(self.data.isf_1[0]).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.isf_create_time_dict["S2"] == data.get("create_time")

        _, file_name = get_cgc_file_name(*self.data.cgc_1[: -1])
        flag, data = CGCInfo(self.data.cgc_1[0]).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.cgc_create_time_dict["S1"] == data.get("create_time")

        # check answer
        for nb in self.isf_num_list:
            isf_tuple = eval("self.data.isf_{}".format(nb))
            isf_param = isf_tuple[: -1]
            isf_answer = isf_tuple[-1]
            if isf_param[0] > end_sn or isf_param[0] < start_sn:
                continue
            _, file_name = get_isf_file_name(*isf_param)
            flag, data = ISFInfo(isf_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data, dict), "self.data.isf_{}, ={} with {}".format(nb, isf_param, isf_answer)
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_isf_file_name(*isf_param, is_full_path=True)
            _, full_finish_file_name = get_isf_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, isf = load_isf(*isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert isf["rows"] == isf_answer["rows"]
            assert isf["cols"] == isf_answer["cols"]
            assert isf["isf"] == isf_answer["isf"], "self.data.isf_{}, ={} with {}".format(nb, isf_param, isf_answer)

            # 以下都是和isf_x σ, μ对称isf的验证
            symmetry_isf_tuple = self.helper.calc_symmetry_isf(*isf_tuple)
            symmetry_isf_param = symmetry_isf_tuple[: -1]
            symmetry_isf_answer = symmetry_isf_tuple[-1]
            flag, symmetry_isf = load_isf(*symmetry_isf_param, is_flag_true_if_not_s_n=True)
            assert flag
            assert symmetry_isf["rows"] == symmetry_isf_answer["rows"]
            assert symmetry_isf["cols"] == symmetry_isf_answer["cols"]
            assert symmetry_isf["isf"] == symmetry_isf_answer["isf"], \
                "self.data.symmetry_isf_{}, ={} with {}".format(nb, symmetry_isf_param, symmetry_isf_answer)

        for nb in self.cgc_num_list:
            cgc_tuple = eval("self.data.cgc_{}".format(nb))
            cgc_param = cgc_tuple[: -1]
            cgc_answer = cgc_tuple[-1]
            if cgc_param[0] > end_sn or cgc_param[0] < start_sn:
                continue
            _, file_name = get_cgc_file_name(*cgc_param)
            flag, data = CGCInfo(cgc_param[0]).query_by_file_name(file_name)
            assert flag
            assert isinstance(data, dict), "self.data.cgc_{}, ={} with {}".format(nb, cgc_param, cgc_answer)
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_cgc_file_name(*cgc_param, is_full_path=True)
            _, full_finish_file_name = get_cgc_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cgc = load_cgc(*cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            cgc_answer_n = cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in cgc_answer.values()) == cgc_answer_n
            assert cgc["N"] == 1
            for cgc_k, cgc_v in cgc_answer.items():
                assert cgc[cgc_k] == Ra(cgc_v) / cgc_answer_n, \
                    "self.data.cgc_{}, ={} with {} with key={} and cgc={}".format(nb, cgc_param, cgc_answer, cgc_k, cgc)
            cgc_answer["N"] = cgc_answer_n

            # 以下都是和cgc_x σ, μ对称cgc的验证
            symmetry_cgc_tuple = self.helper.calc_symmetry_cgc(*cgc_tuple)
            symmetry_cgc_param = symmetry_cgc_tuple[: -1]
            symmetry_cgc_answer = symmetry_cgc_tuple[-1]
            flag, symmetry_cgc = load_cgc(*symmetry_cgc_param, is_flag_true_if_not_s_n=True)
            assert flag
            symmetry_cgc_answer_n = symmetry_cgc_answer.pop("N")
            assert sum(abs(cgc_v) for cgc_v in symmetry_cgc_answer.values()) == symmetry_cgc_answer_n
            assert 1 - isf_0_error_value < abs(symmetry_cgc["N"]) < 1 + isf_0_error_value
            for cgc_k, cgc_v in symmetry_cgc_answer.items():
                assert symmetry_cgc[cgc_k] == Ra(cgc_v) / symmetry_cgc_answer_n, \
                    "symmetry_cgc_{}, ={} with {} with key={} and cgc={}".format(
                        nb, symmetry_cgc_param, symmetry_cgc_answer, cgc_k, cgc)
            symmetry_cgc_answer["N"] = symmetry_cgc_answer_n

        # check finish s_n
        flag, isf_finish_s_n = get_isf_finish_s_n()
        assert flag
        assert isf_finish_s_n == end_sn
        flag, cgc_finish_s_n = get_cgc_finish_s_n()
        assert flag
        assert cgc_finish_s_n == end_sn

        # history_times
        _, isf_finish_file_name = get_isf_finish_s_n_name()
        flag, data = ISFInfo(isf_finish_s_n).query_by_file_name(isf_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 20
        flag, data_txt = ISFInfo(isf_finish_s_n).query_txt_by_file_name(isf_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(2, isf_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

        _, cgc_finish_file_name = get_cgc_finish_s_n_name()
        flag, data = CGCInfo(cgc_finish_s_n).query_by_file_name(cgc_finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 100
        flag, data_txt = CGCInfo(cgc_finish_s_n).query_txt_by_file_name(cgc_finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(1, cgc_finish_s_n + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)

    '''

    def test_data_helper_init(self):
        # 检查非法范围报错
        flag = True
        try:
            data_s_0 = DataHelper(0)
            flag = False
        except Exception as e:
            assert str(e) == "s_k={} must be int and >= 1".format(0)
        assert flag

        # 检查最小初始情况
        data_sn = DataHelper(1)
        assert data_sn.s_n == 1
        assert data_sn.yd_list == [[1]]
        assert data_sn.bl_yd_list_dict == {tuple([1]): [[]]}
        assert data_sn.yt_num_dict == {tuple([1]): 1}
        assert data_sn.eigenvalue_list == [1]

        # 检查传递和更新
        data_st = data_sn
        data_sn = DataHelper(2)
        assert data_st.s_n == 1
        assert data_st.yd_list == [[1]]
        assert data_st.bl_yd_list_dict == {tuple([1]): [[]]}
        assert data_st.yt_num_dict == {tuple([1]): 1}
        assert data_st.eigenvalue_list == [1]

        assert data_sn.s_n == 2
        assert data_sn.yd_list == [[2], [1, 1]]
        assert data_sn.bl_yd_list_dict == {tuple([2]): [[1]], tuple([1, 1]): [[1]]}
        assert data_sn.yt_num_dict == {tuple([2]): 1, tuple([1, 1]): 1}
        assert data_sn.eigenvalue_list == [1, -1]

    def test_σ_μ_data_helper_init(self):
        # 检查最小初始情况
        data_st = DataHelper(1)
        data_sn = DataHelper(2)
        data_σ_μ = ΣMDataHelper(2, [2], [1, 1], data_sn, data_st)
        assert data_σ_μ.s_n == 2
        assert data_σ_μ.s_t == 1
        assert data_σ_μ.σ == [2]
        assert data_σ_μ.μ == [1, 1]
        assert data_σ_μ.data_sn_cls == data_sn
        assert data_σ_μ.data_st_cls == data_st
        assert (data_σ_μ.cg_series_list == np.array([0, 1], dtype=int)).all()
        assert data_σ_μ.bl_yds_of_σ == [[1]]
        assert data_σ_μ.bl_yds_of_μ == [[1]]
        assert data_σ_μ.cg_series_st_list_dict == {(tuple([1]), tuple([1])): [1]}
        assert data_σ_μ.in_matrix_σ_dict == {(1, 2): sp.Matrix([[1]])}
        assert data_σ_μ.in_matrix_μ_dict == {(1, 2): sp.Matrix([[-1]])}

        # 再来一个
        data_st = DataHelper(3)
        data_sn = DataHelper(4)
        data_σ_μ = ΣMDataHelper(4, [3, 1], [2, 2], data_sn, data_st)
        assert data_σ_μ.s_n == 4
        assert data_σ_μ.s_t == 3
        assert data_σ_μ.σ == [3, 1]
        assert data_σ_μ.μ == [2, 2]
        assert data_σ_μ.data_sn_cls == data_sn
        assert data_σ_μ.data_st_cls == data_st
        assert (data_σ_μ.cg_series_list == np.array([0, 1, 0, 1, 0], dtype=int)).all()
        assert data_σ_μ.bl_yds_of_σ == [[3], [2, 1]]
        assert data_σ_μ.bl_yds_of_μ == [[2, 1]]
        cg_series_st_list_dict_answer = {(tuple([3]), tuple([2, 1])): np.array([0, 1, 0], dtype=int),
                                         (tuple([2, 1]), tuple([2, 1])): np.array([1, 1, 1], dtype=int)}
        for k, v in data_σ_μ.cg_series_st_list_dict.items():
            assert (v == cg_series_st_list_dict_answer.get(k)).all()

        in_matrix_σ_dict_answer = {(1, 4): sp.Matrix([[-Ra(1)/3, -sqrt(2)/3, -sqrt(6)/3],
                                                      [-sqrt(2)/3, Ra(5)/6, -Ra(1)/sqrt(12)],
                                                      [-sqrt(6)/3, -Ra(1)/sqrt(12), Ra(1)/2]]),
                                   (2, 4): sp.Matrix([[-Ra(1)/3, -sqrt(2)/3, sqrt(6)/3],
                                                      [-sqrt(2)/3, Ra(5)/6, Ra(1)/sqrt(12)],
                                                      [sqrt(6)/3, Ra(1)/sqrt(12), Ra(1)/2]]),
                                   (3, 4): sp.Matrix([[-Ra(1)/3, sqrt(8)/3, 0],
                                                      [sqrt(8)/3, Ra(1)/3, 0],
                                                      [0, 0, 1]])}
        for k, v in data_σ_μ.in_matrix_σ_dict.items():
            assert v == in_matrix_σ_dict_answer.get(k)

        in_matrix_μ_dict_answer = {(1, 4): sp.Matrix([[-Ra(1)/2, sqrt(3)/2],
                                                      [sqrt(3)/2, Ra(1)/2]]),
                                   (2, 4): sp.Matrix([[-Ra(1)/2, -sqrt(3)/2],
                                                      [-sqrt(3)/2, Ra(1)/2]]),
                                   (3, 4): sp.Matrix([[1, 0],
                                                      [0, -1]])}
        for k, v in data_σ_μ.in_matrix_μ_dict.items():
            assert v == in_matrix_μ_dict_answer.get(k)

    def test_calc_helper(self):
        # 检查非法范围报错
        calc_helper = CalcHelper()
        flag = True
        try:
            calc_helper.enable_now_s_n(0)
            flag = False
        except Exception as e:
            assert str(e) == "s_n={} must be int and >= {}".format(0, 1)
        assert flag

        # 检查enable_now_s_n
        calc_helper = CalcHelper()
        calc_helper.enable_now_s_n(3)
        assert calc_helper.s_n == 3
        assert calc_helper.s_t == 2

        # 检查_calc_m_with_m_st
        assert 2 == calc_helper._calc_m_with_m_st([2, 1], 2, [[2, 1]], {tuple([2, 1]): 3})

        data_st = DataHelper(4)
        data_sn = DataHelper(5)
        bl_yd_list = data_sn.bl_yd_list_dict[tuple([3, 2])]
        assert 4 == calc_helper._calc_m_with_m_st([2, 2], 1, bl_yd_list, data_st.yt_num_dict)

        if self.test_sn >= 6:
            data_st = data_sn
            data_sn = DataHelper(6)
            bl_yd_list = data_sn.bl_yd_list_dict[tuple([3, 2, 1])]
            assert 14 == calc_helper._calc_m_with_m_st([2, 2, 1], 3, bl_yd_list, data_st.yt_num_dict)
            assert 0 == calc_helper._calc_m_with_m_st([3, 2], 0, bl_yd_list, data_st.yt_num_dict)
            assert 5 == calc_helper._calc_m_with_m_st([3, 1, 1], 0, bl_yd_list, data_st.yt_num_dict)
            assert 11 == calc_helper._calc_m_with_m_st([2, 2, 1], 0, bl_yd_list, data_st.yt_num_dict)

            bl_yd_list = data_sn.bl_yd_list_dict[tuple([2, 2, 1, 1])]
            assert 5 == calc_helper._calc_m_with_m_st([2, 1, 1, 1], 0, bl_yd_list, data_st.yt_num_dict)

    def test_isf_helper_init(self):
        isf_func = ISFHelper()
        isf_func.enable_now_s_n(2)
        assert callable(isf_func.calc_isf_dict)
        # assert callable(isf_func._calc_isf_phase)
        # assert callable(isf_func._calc_before_isf_ref)
        # assert callable(isf_func._calc_ref_isf_tmp)
        assert callable(isf_func._calc_isf_β_and_phase_list)
        assert callable(isf_func._calc_isf_fbl_another_way)
        assert callable(isf_func._get_first_no_0_number_from_vector)
        assert callable(isf_func._calc_orthogonalization_vector)
        assert callable(isf_func._calc_schmidt_orthogonalization_tricky)
        # assert callable(isf_func._calc_schmidt_orthogonalization_eigenvectors_and_β_list)
        assert callable(isf_func._calc_ν_by_λ_and_bl)
        assert callable(isf_func.calc_row_indexes_tmp)
        assert callable(isf_func._load_cgc_with_m1_by_input_json)
        assert callable(isf_func._load_cgc_with_m1_lru)
        assert callable(isf_func._cgc_st_2_cgc_m_dict)
        assert callable(isf_func._calc_isf_matrix_element)
        assert callable(isf_func._calc_isf_matrix)

    def test_isf_helper_calc_row_indexes_tmp(self):
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, β'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        isf_func = ISFHelper()

        # 小
        data_st = DataHelper(1)
        data_sn = DataHelper(2)
        data_σ_μ = ΣMDataHelper(2, [2], [1, 1], data_sn, data_st)
        answer = [([1], [1], None)]
        row_index_tmp_list = isf_func.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                           data_σ_μ.cg_series_st_list_dict, [1], data_st.yd_list)
        assert answer == row_index_tmp_list

        # 大
        data_st = DataHelper(4)
        data_sn = DataHelper(5)
        data_σ_μ = ΣMDataHelper(5, [2, 1, 1, 1], [2, 1, 1, 1], data_sn, data_st)
        answer = [([2, 1, 1], [2, 1, 1], None), ([1, 1, 1, 1], [1, 1, 1, 1], None)]
        row_index_tmp_list = isf_func.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                           data_σ_μ.cg_series_st_list_dict,
                                                           [4], data_st.yd_list)
        assert answer == row_index_tmp_list

        answer = []
        row_index_tmp_list = isf_func.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                           data_σ_μ.cg_series_st_list_dict,
                                                           [1, 1, 1, 1], data_st.yd_list)
        assert answer == row_index_tmp_list

        # 有多重
        if self.test_sn >= 6:
            data_st = DataHelper(5)
            data_sn = DataHelper(6)
            data_σ_μ = ΣMDataHelper(6, [3, 2, 1], [3, 2, 1], data_sn, data_st)
            answer = [([3, 2], [3, 2], None), ([3, 2], [3, 1, 1], None), ([3, 2], [2, 2, 1], None),
                      ([3, 1, 1], [3, 2], None), ([3, 1, 1], [3, 1, 1], 1),
                      ([3, 1, 1], [3, 1, 1], 2), ([3, 1, 1], [2, 2, 1], None),
                      ([2, 2, 1], [3, 2], None), ([2, 2, 1], [3, 1, 1], None), ([2, 2, 1], [2, 2, 1], None)]
            cols_index = [([4, 2], 1), ([4, 2], 2), ([4, 2], 3),
                          ([3, 3], 1), ([3, 3], 2),
                          ([3, 2, 1], 1), ([3, 2, 1], 2), ([3, 2, 1], 3), ([3, 2, 1], 4), ([3, 2, 1], 5)]
            row_index_tmp_list = isf_func.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                               data_σ_μ.cg_series_st_list_dict,
                                                               [3, 2], data_st.yd_list)
            assert answer == row_index_tmp_list

    def test_isf_helper__calc_ν_by_λ_and_bl(self):
        isf_func = ISFHelper()

        # 小值
        data_sn = DataHelper(1)
        assert [1] == isf_func._calc_ν_by_λ_and_bl(1, data_sn.eigenvalue_list, data_sn.yd_list, None, None)
        data_sn = DataHelper(2)
        assert [2] == isf_func._calc_ν_by_λ_and_bl(1, data_sn.eigenvalue_list, data_sn.yd_list, None, None)
        assert [1, 1] == isf_func._calc_ν_by_λ_and_bl(-1, data_sn.eigenvalue_list, data_sn.yd_list, None, None)

        if self.test_sn >= 6:
            # 无偶然简并
            data_sn = DataHelper(6)
            assert [6] == isf_func._calc_ν_by_λ_and_bl(15, data_sn.eigenvalue_list, data_sn.yd_list, None, None)
            assert [3, 2, 1] == isf_func._calc_ν_by_λ_and_bl(0, data_sn.eigenvalue_list, data_sn.yd_list, None, None)
            assert [2, 1, 1, 1, 1] == isf_func._calc_ν_by_λ_and_bl(-9, data_sn.eigenvalue_list, data_sn.yd_list, None, None)
            # 有偶然简并
            assert [3, 3] == isf_func._calc_ν_by_λ_and_bl(3, data_sn.eigenvalue_list, data_sn.yd_list,
                                                          data_sn.bl_yd_list_dict, [3, 2])
            assert [4, 1, 1] == isf_func._calc_ν_by_λ_and_bl(3, data_sn.eigenvalue_list, data_sn.yd_list,
                                                             data_sn.bl_yd_list_dict, [4, 1])
            assert [4, 1, 1] == isf_func._calc_ν_by_λ_and_bl(3, data_sn.eigenvalue_list, data_sn.yd_list,
                                                             data_sn.bl_yd_list_dict, [3, 1, 1])

    def test_isf_helper__calc_orthogonalization_vector(self):
        isf_func = ISFHelper()
        vector_list = [sp.Matrix([sqrt(3)]),
                       sp.Matrix([-sqrt(3)]),
                       sp.Matrix([1, 1]),
                       sp.Matrix([-1, 1]),
                       sp.Matrix([-1, 1]),
                       sp.Matrix([Ra(1)/sqrt(2), Ra(1)/sqrt(3), Ra(1)/sqrt(5), 0]),
                       sp.Matrix([sqrt(7), sqrt(7), sqrt(7)])]
        answer_list = [sp.Matrix([1]),
                       sp.Matrix([-1]),
                       sp.Matrix([sqrt(2)/2, sqrt(2)/2]),
                       sp.Matrix([-sqrt(2)/2, sqrt(2)/2]),
                       sp.Matrix([-sqrt(2)/2, -sqrt(2)/-2]),
                       sp.Matrix([sqrt(30)/sqrt(31*2), sqrt(30)/sqrt(31*3), sqrt(30)/sqrt(31*5), 0]),
                       sp.Matrix([sqrt(3)/3, sqrt(3)/3, sqrt(3)/3])]
        for v, a in zip(vector_list, answer_list):
            assert a == isf_func._calc_orthogonalization_vector(v)

    def test_isf_helper__calc_schmidt_orthogonalization_tricky(self):
        pass

    def test_isf_helper__get_first_no_0_number_from_vector(self):
        isf_func = ISFHelper()
        assert (None, None) == isf_func._get_first_no_0_number_from_vector([0, 0, 0, 0])
        # assert isf_func._get_first_no_0_number_from_vector([0.0000001, 0.0000001, 0.0000001]) is None
        # assert isf_func._get_first_no_0_number_from_vector([0.0000009, -0.0000009, 0.0000009]) is None
        assert (1, 1) == isf_func._get_first_no_0_number_from_vector([0, 1, 0, 0])
        assert (2, -3) == isf_func._get_first_no_0_number_from_vector([0, 0, -3, 0])
        assert (1, 1) == isf_func._get_first_no_0_number_from_vector([0, 1, 1, -3])
        assert (0, 0.0000009) == isf_func._get_first_no_0_number_from_vector(sp.Matrix([0.0000009, -0.10010, 0.50009]))

    def test_isf_helper__calc_before_isf_reference(self):
        pass

    def test_isf_helper__calc_isf_phase(self):
        pass

    def test_isf_helper_calc_isf_dict(self):
        pass

    def test_isf_helper__load_cgc_with_m1_by_input_json(self):
        pass

    def test_isf_helper__load_cgc_with_m1_lru(self):
        pass

    def test_isf_helper__cgc_st_2_cgc_m_dict(self):
        pass

    def test_isf_helper__calc_isf_matrix_element(self):
        pass

    def test_isf_helper__calc_isf_matrix(self):
        pass

    def test_cgc_helper_init(self):
        cgc_func = CGCHelper()
        cgc_func.enable_now_s_n(1)
        assert callable(cgc_func.calc_cgc_dict_part_and_save_by_isf)
        
    '''