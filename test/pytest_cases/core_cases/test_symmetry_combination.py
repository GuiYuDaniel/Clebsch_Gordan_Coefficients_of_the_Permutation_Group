# -*- coding:utf8 -*-
"""
测试
core/symmetry_combination.py
下所有功能是否正确执行
"""


import os
import pytest
import time
from conf.cgc_config import default_s_n, cgc_db_name
from core.cgc_utils.cgc_db_typing import SYMInfo
from core.cgc_utils.cgc_local_db import get_meta_σμν_file_name, get_sym_σμν_file_name
from core.cgc_utils.cgc_local_db import get_symmetry_combination_finish_s_n_name
from core.young_diagrams import create_young_diagrams, load_young_diagrams
from core.characters_and_gi import create_characters_and_gi
from core.cg_series import create_cg_series
from core.eigenvalues import create_eigenvalues
from core.symmetry_combination import load_meta_σμν, load_sym_σμν
from core.symmetry_combination import get_symmetry_combination_finish_s_n, create_symmetry_combination
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestSYM(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_db_name, extension_name=".test_young_diagrams_protected")
        self.protector.protector_setup()

        # 准备前文
        s_n = 9
        self.test_sn = s_n
        flag, msg = create_young_diagrams(s_n)
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

        self.meta_σμν_tuple_list_s_1 = [([1], [1], [1])]
        self.meta_σμν_tuple_list_s_3 = [([3], [3], [3]), ([2, 1], [2, 1], [3]), ([2, 1], [2, 1], [2, 1])]

        self.sym_σμν_dict_s_1 = {('[1]', '[1]', '[1]'): ['σμν', 'σ~μ~ν', 'σ~μν~', 'σμ~ν~',
                                                         'μσν', 'μ~σ~ν', 'μ~σν~', 'μσ~ν~',
                                                         'νμσ', 'ν~μ~σ', 'ν~μσ~', 'νμ~σ~',
                                                         'σνμ', 'σ~ν~μ', 'σ~νμ~', 'σν~μ~',
                                                         'νσμ', 'ν~σ~μ', 'ν~σμ~', 'νσ~μ~',
                                                         'μνσ', 'μ~ν~σ', 'μ~νσ~', 'μν~σ~']}
        self.sym_σμν_dict_s_3 = {('[3]', '[3]', '[3]'): ['σμν', 'μσν', 'νμσ', 'σνμ', 'νσμ', 'μνσ'],
                                 ('[1, 1, 1]', '[1, 1, 1]', '[3]'): ['σ~μ~ν', 'μ~σ~ν', 'ν~μ~σ',
                                                                     'σ~ν~μ', 'ν~σ~μ', 'μ~ν~σ'],
                                 ('[1, 1, 1]', '[3]', '[1, 1, 1]'): ['σ~μν~', 'μ~σν~', 'ν~μσ~',
                                                                     'σ~νμ~', 'ν~σμ~', 'μ~νσ~'],
                                 ('[3]', '[1, 1, 1]', '[1, 1, 1]'): ['σμ~ν~', 'μσ~ν~', 'νμ~σ~',
                                                                     'σν~μ~', 'νσ~μ~', 'μν~σ~']}

        self.default_s_n = default_s_n
        _, self.s_1_meta_file_name = get_meta_σμν_file_name(1)
        _, self.s_1_meta_full_file_name = get_meta_σμν_file_name(1, is_full_path=True)
        _, self.s_3_meta_file_name = get_meta_σμν_file_name(3)
        _, self.s_3_meta_full_file_name = get_meta_σμν_file_name(3, is_full_path=True)
        _, self.s_5_meta_file_name = get_meta_σμν_file_name(5)
        _, self.s_5_meta_full_file_name = get_meta_σμν_file_name(5, is_full_path=True)
        _, self.s_1_sym_file_name = get_sym_σμν_file_name(1, [1], [1], [1])
        _, self.s_1_sym_full_file_name = get_sym_σμν_file_name(1, [1], [1], [1], is_full_path=True)
        _, self.s_3_sym_file_name = get_sym_σμν_file_name(3, [3], [3], [3])
        _, self.s_3_sym_full_file_name = get_sym_σμν_file_name(3, [3], [3], [3], is_full_path=True)
        _, self.s_n_finish_file_name = get_symmetry_combination_finish_s_n_name()
        _, self.s_n_finish_meta_full_file_name = get_symmetry_combination_finish_s_n_name(is_full_path=True)

        self.create_time_dict = {}  # 用于检查计算好的部分不会重复计算

    def teardown_class(self):
        # self.protector.protector_teardown()
        pass

        # start with 0xx tests need test by order

    def test_001_create_symmetry_combination_s_n_1(self):  # calc and save, return True, None
        """for s_n=1, there is no finish db"""
        # check with no db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.s_1_meta_full_file_name + ex)
            assert not os.path.exists(self.s_n_finish_meta_full_file_name + ex)
        flag, meta_σμν_tuple_list = load_meta_σμν(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert meta_σμν_tuple_list is False
        flag, meta_σμν_tuple_list = load_meta_σμν(1, is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(meta_σμν_tuple_list, str)
        flag, finish_s_n = get_symmetry_combination_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # check create_young_diagrams_s_n_1
        flag, msg = create_symmetry_combination(1)
        assert flag
        assert msg == 1

        # check answer
        flag, data = SYMInfo(1).query_by_file_name(self.s_1_meta_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        self.create_time_dict["S1"] = data.get("create_time")

        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_1_meta_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_meta_full_file_name + ex)
        flag, meta_σμν_tuple_list = load_meta_σμν(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert meta_σμν_tuple_list == self.meta_σμν_tuple_list_s_1

        flag, data = SYMInfo(1).query_by_file_name(self.s_1_sym_file_name)
        assert flag
        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_1_sym_full_file_name + ex)
        flag, sym_σμν_dict = load_sym_σμν(1, [1], [1], [1], is_flag_true_if_not_s_n=True)
        assert flag
        assert sym_σμν_dict == self.sym_σμν_dict_s_1

        # check finish s_n
        flag, finish_s_n = get_symmetry_combination_finish_s_n()
        assert flag
        assert finish_s_n == 1

        # history_times
        _, finish_file_name = get_symmetry_combination_finish_s_n_name()
        flag, data = SYMInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        assert isinstance(data.get("flags").get("history_times").get("S1"), int)
        assert 0 <= data.get("flags").get("history_times").get("S1") <= 1
        flag, data_txt = SYMInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        assert isinstance(data.get("history_times").get("S{}".format(1)), int)

    def test_002_create_symmetry_combination_s_n_2_to_5(self):
        # check create_young_diagrams_s_n 2 to 5
        flag, msg = create_symmetry_combination(5)
        assert flag
        assert msg == 5

        # check answer
        flag, data = SYMInfo(3).query_by_file_name(self.s_3_meta_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        self.create_time_dict["S3"] = data.get("create_time")
        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_3_meta_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_meta_full_file_name + ex)
        flag, meta_σμν_tuple_list = load_meta_σμν(3, is_flag_true_if_not_s_n=True)
        assert flag
        assert meta_σμν_tuple_list == self.meta_σμν_tuple_list_s_3

        flag, data = SYMInfo(3).query_by_file_name(self.s_3_sym_file_name)
        assert flag
        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_3_sym_full_file_name + ex)
        flag, sym_σμν_dict = load_sym_σμν(3, [3], [3], [3], is_flag_true_if_not_s_n=True)
        assert flag
        assert sym_σμν_dict == self.sym_σμν_dict_s_3

        flag, data = SYMInfo(1).query_by_file_name(self.s_1_meta_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get("S1") == data.get("create_time")

        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_5_meta_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_meta_full_file_name + ex)
        flag, finish_s_n = get_symmetry_combination_finish_s_n()
        assert flag
        assert finish_s_n == 5
        _, finish_file_name = get_symmetry_combination_finish_s_n_name()
        flag, data = SYMInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2, 3, 4, 5]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
        for i in [6, 7]:
            assert data.get("flags").get("history_times").get("S{}".format(i)) is None
        flag, data_txt = SYMInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4, 5]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        for i in [6, 7]:
            assert data.get("history_times").get("S{}".format(i)) is None
