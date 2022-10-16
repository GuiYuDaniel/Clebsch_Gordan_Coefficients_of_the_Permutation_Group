# -*- coding:utf8 -*-
"""
测试
core/characters_and_gi.py
下所有功能是否正确执行
"""


# 1，用历史数据检查；
# 2，用正交归一性检查；
# 3，用分支律生成检查；  # TODO
# 4，用（-1^n检查）；
# 5，用len检查；
# 6，用生成的杨图检查
# 7，用sum(gi) == n!检查
import os
import time
import math
import numpy as np
from conf.cgc_config import cgc_rst_folder
from core.branching_laws import create_branching_laws
from core.young_diagrams import create_young_diagrams, load_young_diagrams, calc_young_diagram_tilde
from core.characters_and_gi import CharacterData, GiData
from core.characters_and_gi import load_characters_and_gi, get_characters_and_gi_finish_s_n
from core.characters_and_gi import create_characters_and_gi, calc_single_characters_and_gi
from core.cgc_utils.cgc_db_typing import CharacterAndGiInfo
from core.cgc_utils.cgc_local_db import get_characters_and_gi_file_name, get_characters_and_gi_finish_s_n_name
from test.fake.test_results.fake_cgc.old_characters_and_gi_data import *
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestCharacterData(object):
    """关于录入数据，必须尽所有可能去检查！"""

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_characters_protected")
        self.protector.protector_setup()
        self.new_data_cls = CharacterData()
        self.new_gi_cls = GiData()

        # 准备前文
        flag, msg = create_young_diagrams(9)
        assert flag
        assert msg == 9
        flag, msg = create_branching_laws(9)
        assert flag
        assert msg == 9

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    def test_by_old_data(self):
        # test from S3 to S7 by old data
        for s_i in range(3, 7+1):
            old_character_dict = eval("matrix_{}".format(s_i))
            new_character_dict = eval("self.new_data_cls.character_{}".format(s_i))
            flag, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
            assert flag
            yd_len = len(yd_list)
            assert len(old_character_dict) == yd_len
            assert len(new_character_dict) == yd_len
            old_ch_list = [old_character_dict[tuple(i)] for i in yd_list]
            new_ch_list = [new_character_dict[tuple(i)] for i in yd_list]
            old_ch_matrix = np.array(old_ch_list, dtype=int)
            new_ch_matrix = np.array(new_ch_list, dtype=int)
            assert old_ch_matrix.size == new_ch_matrix.size
            old_ch_t_list = old_ch_matrix.T.tolist()
            new_ch_t_list = new_ch_matrix.T.tolist()
            for single_list in new_ch_t_list:
                assert single_list in old_ch_t_list

        # test gi
        for s_i in range(3, 7 + 1):
            old_gi_list = eval("gi_{}".format(s_i))
            new_gi_list = eval("self.new_gi_cls.gi_{}".format(s_i))
            assert len(old_gi_list) == len(new_gi_list)
            for g_i in new_gi_list:
                assert g_i in old_gi_list

    def test_by_others(self):
        # test from all character which from S1 to S9 by property
        all_ch_matrix_list = []
        for s_i in range(1, 9 + 1):
            # logger.warning("@@@@ s_i={}".format(s_i))
            new_character_dict = eval("self.new_data_cls.character_{}".format(s_i))
            new_gi_list = eval("self.new_gi_cls.gi_{}".format(s_i))
            flag, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
            assert flag
            matrix_div = len(yd_list)
            assert len(new_character_dict) == matrix_div  # 检查字典的维度
            new_ch_list = [new_character_dict[tuple(i)] for i in yd_list]  # 用生成的杨图检查
            # new_ch_len_list = [len(i) for i in new_ch_list]
            # logger.warning("@@@@ new_ch_list={} with new_ch_len_list={}".format(new_ch_list, new_ch_len_list))
            # if s_i == 9:
            #     logger.warning("$$$$ 7={} len(7)={}".format(new_ch_list[7], len(new_ch_list[7])))
            new_ch_matrix = np.array(new_ch_list, dtype=int)  # 这里可以检查列的维度
            new_gi_array = np.array(new_gi_list, dtype=int)  # TODO 小心！！！dtype会把小数取整！
            # logger.warning("@@@@ new_ch_matrix={}".format(new_ch_matrix))
            last_array = new_ch_matrix[-1, :]
            middle_s_i = s_i / 2

            # 用（-1^n检查）
            for i in range(matrix_div):
                yd = yd_list[i]
                if yd[0] <= middle_s_i:
                    i_array = new_ch_matrix[i, :]
                    _, yd_tilde = calc_young_diagram_tilde(yd)
                    if yd_tilde == yd:
                        continue
                    # 剩下的都是 首行 <= Sn/2 又不自共轭的，它们都可以用（-1^n检查）
                    yd_tilde_index = yd_list.index(yd_tilde)
                    i_tilde_array = new_ch_matrix[yd_tilde_index, :]
                    assert (i_tilde_array == i_array * last_array).all(), \
                        "with yd={}, yd_tilde={}".format(yd, yd_tilde)

            # 用sum(gi) == n!检查
            big_g = sum(new_gi_list)
            assert big_g == math.factorial(s_i)

            # 用正交归一性检查
            # logger.warning("@@@@ new_ch_matrix={}".format(new_ch_matrix))
            for i in range(matrix_div):
                # yd = yd_list[i]
                # _, yd_tilde = calc_young_diagram_tilde(yd)
                # yd_tilde_index = yd_list.index(yd_tilde)
                for j in range(i, matrix_div):

                    # 列正交
                    if i == j:
                        continue
                    assert sum(new_ch_matrix[:, i] * new_ch_matrix[:, j]) == 0, \
                        "{} * {} should be 0".format(new_ch_matrix[:, i], new_ch_matrix[:, j])

                    # 行正交归G
                    sum_i_dot_j_dot_g = sum(new_ch_matrix[i, :] * new_ch_matrix[j, :] * new_gi_array)
                    if j == i:
                        assert sum_i_dot_j_dot_g == big_g, \
                            "{} * {} * {} should be big_g".format(new_ch_matrix[i, :],
                                                                  new_ch_matrix[j, :], new_gi_array)
                    else:
                        assert sum_i_dot_j_dot_g == 0, \
                            "{} * {} * {} should be 0".format(new_ch_matrix[i, :], new_ch_matrix[j, :], new_gi_array)

            all_ch_matrix_list.append(new_ch_matrix)


class TestCharacters(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_characters_protected")
        self.protector.protector_setup()

        # 准备前文
        flag, msg = create_young_diagrams(6)
        assert flag
        assert msg == 6

        self.characters_and_gi_s_1 = {"character": np.array([[1]]),
                                      "gi": np.array([1])}
        self.characters_and_gi_s_2 = {"character": np.array([[1, 1], [1, -1]]),
                                      "gi": np.array([1, 1])}
        self.characters_and_gi_s_3 = {"character": np.array([[1, 1, 1], [2, 0, -1], [1, -1, 1]]),
                                      "gi": np.array([1, 3, 2])}
        self.characters_and_gi_s_4 = {"character": np.array([[1, 1, 1, 1, 1], [3, 1, 0, -1, -1], [2, 0, -1, 0, 2],
                                                             [3, -1, 0, 1, -1], [1, -1, 1, -1, 1]]),
                                      "gi": np.array([1, 6, 8, 6, 3])}
        self.characters_and_gi_s_5 = {"character": np.array([[1, 1, 1, 1, 1, 1, 1], [4, 2, 1, 0, 0, -1, -1],
                                                             [5, 1, -1, -1, 1, 1, 0], [6, 0, 0, 0, -2, 0, 1],
                                                             [5, -1, -1, 1, 1, -1, 0], [4, -2, 1, 0, 0, 1, -1],
                                                             [1, -1, 1, -1, 1, -1, 1]]),
                                      "gi": np.array([1, 10, 20, 30, 15, 20, 24])}

        _, self.s_1_file_name = get_characters_and_gi_file_name(1)
        _, self.s_1_full_file_name = get_characters_and_gi_file_name(1, is_full_path=True)
        _, self.s_n_finish_file_name = get_characters_and_gi_finish_s_n_name()
        _, self.s_n_finish_full_file_name = get_characters_and_gi_finish_s_n_name(is_full_path=True)

        _, self.s_4_full_file_name = get_characters_and_gi_file_name(4, is_full_path=True)
        self.create_time_dict = {}  # 用于检查计算好的部分不会重复计算

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    def test_001_create_characters_and_gi(self):
        """for s_n=1, there is no finish db"""
        # check with no db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.s_1_full_file_name + ex)
            assert not os.path.exists(self.s_n_finish_full_file_name + ex)
        flag, characters_and_gi = load_characters_and_gi(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert characters_and_gi is False
        flag, characters_and_gi = load_characters_and_gi(1, is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(characters_and_gi, str)
        flag, finish_s_n = get_characters_and_gi_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # check create_characters_and_gi_s_n_1
        flag, msg = create_characters_and_gi(1)
        assert flag
        assert msg == 1

        # check answer
        flag, data = CharacterAndGiInfo(1).query_by_file_name(self.s_1_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        self.create_time_dict["S1"] = data.get("create_time")
        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_1_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_full_file_name + ex)
        flag, characters_and_gi = load_characters_and_gi(1, is_flag_true_if_not_s_n=True)
        assert flag
        assert characters_and_gi == self.characters_and_gi_s_1

        # check finish s_n
        flag, finish_s_n = get_characters_and_gi_finish_s_n()
        assert flag
        assert finish_s_n == 1

        # history_times
        _, finish_file_name = get_characters_and_gi_finish_s_n_name()
        flag, data = CharacterAndGiInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        assert isinstance(data.get("flags").get("history_times").get("S1"), int)
        assert 0 <= data.get("flags").get("history_times").get("S1") <= 1
        flag, data_txt = CharacterAndGiInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        assert isinstance(data.get("history_times").get("S{}".format(1)), int)

    def test_002_create_characters_and_gi_s_n_2_to_4(self):
        # check create_characters_and_gi from S2 to S4
        flag, msg = create_characters_and_gi(4)
        assert flag
        assert msg == 4

        flag, data = CharacterAndGiInfo(1).query_by_file_name(self.s_1_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get("S1") == data.get("create_time")

        for ex in [".pkl", ".txt"]:
            assert os.path.exists(self.s_4_full_file_name + ex)
            assert os.path.exists(self.s_n_finish_full_file_name + ex)
        for i in [1, 2, 3, 4]:
            flag, characters_and_gi = load_characters_and_gi(i, is_flag_true_if_not_s_n=True)
            assert flag
            assert (characters_and_gi["character"] == eval("self.characters_and_gi_s_{}".format(i))["character"]).all()
            assert (characters_and_gi["gi"] == eval("self.characters_and_gi_s_{}".format(i))["gi"]).all()
            _, file_name = get_characters_and_gi_file_name(i)
            _, data = CharacterAndGiInfo(i).query_by_file_name(file_name)
            _, yd_list_s_i = load_young_diagrams(i, is_flag_true_if_not_s_n=False)
            assert isinstance(data, dict)
            assert data.get("flags").get("young_diagram_index") == yd_list_s_i
        flag, finish_s_n = get_characters_and_gi_finish_s_n()
        assert flag
        assert finish_s_n == 4
        _, finish_file_name = get_characters_and_gi_finish_s_n_name()
        flag, data = CharacterAndGiInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == {}
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("flags").get("history_times").get("S{}".format(i)) is None
        flag, data_txt = CharacterAndGiInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("history_times").get("S{}".format(i)) is None

    def test_calc_single_characters_and_gi(self):
        """1 to 4 is finish, 5 is needing calc"""
        head_s_n = 1
        tail_s_n = 5
        for i in range(head_s_n, tail_s_n + 1):
            flag, characters_and_gi = calc_single_characters_and_gi(i)
            assert flag
            assert (characters_and_gi["character"] == eval("self.characters_and_gi_s_{}".format(i))["character"]).all()
            assert (characters_and_gi["gi"] == eval("self.characters_and_gi_s_{}".format(i))["gi"]).all()
        _, finish_file_name = get_characters_and_gi_finish_s_n_name()
        flag, data = CharacterAndGiInfo(0).query_by_file_name(finish_file_name)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("flags").get("history_times").get("S{}".format(i)) is None
        flag, data_txt = CharacterAndGiInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in [1, 2, 3, 4]:
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        for i in [5, 6, 7]:
            assert data.get("history_times").get("S{}".format(i)) is None
