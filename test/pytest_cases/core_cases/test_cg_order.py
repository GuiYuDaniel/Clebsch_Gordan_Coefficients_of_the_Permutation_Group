# -*- coding:utf8 -*-
"""
测试
core/cg_series.py
下所有功能是否正确执行
"""


import os
import time
import math
import numpy as np
from conf.cgc_config import cgc_rst_folder
from core.young_diagrams import create_young_diagrams
from core.characters_and_gi import create_characters_and_gi
from core.cg_series import create_cg_series, save_cg_series, load_cg_series
from core.cg_series import save_cg_series_finish_s_n, get_cg_series_finish_s_n
from core.cgc_utils.cgc_db_typing import CGSeriesInfo
from core.cgc_utils.cgc_local_db import get_cg_series_file_name, get_cg_series_finish_s_n_name
from db.local_db_protector import DBProtector
from utils.log import get_logger


logger = get_logger(__name__)


class TestCGSeries(object):

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cg_series_protected")
        self.protector.protector_setup()

        # 准备前文
        flag, msg = create_young_diagrams(6)
        assert flag
        assert msg == 6

        flag, msg = create_characters_and_gi(6)
        assert flag
        assert msg == 6

        # by old data
        # 格式(s_n, yd_1, yd_2, cg_series_list)
        self.cgo_case_1 = (1, [1], [1], [1])
        self.cgo_case_2 = (2, [2], [2], [1, 0])
        self.cgo_case_3 = (2, [2], [1, 1], [0, 1])
        self.cgo_case_4 = (2, [1, 1], [2], [0, 1])
        self.cgo_case_5 = (2, [1, 1], [1, 1], [1, 0])
        self.cgo_case_6 = (3, [3], [2, 1], [0, 1, 0])
        self.cgo_case_7 = (4, [1, 1, 1, 1], [2, 1, 1], [0, 1, 0, 0, 0])
        self.cgo_case_8 = (5, [3, 2], [5], [0, 0, 1, 0, 0, 0, 0])
        self.cgo_case_9 = (6, [2, 2, 1, 1], [2, 2, 1, 1], [1, 1, 2, 1, 0, 2, 1, 1, 0, 0, 0])
        self.cgo_case_10 = (7, [4, 3], [5, 2], [0, 1, 1, 1, 1, 2, 0, 1, 1, 1, 0, 0, 0, 0, 0])
        self.cgo_case_num_list = list(range(1, 10 + 1))  # [1, ..., 10]

        # by books
        # 格式(yd_1_tuple, yd_2_tuple, bool_tuple, cg_series_list, s_n)
        self.cgo_book_1 = (([2, 1],), ([2, 1],), (True,), [1, 1, 1], 3)
        self.cgo_book_2 = (([3, 1], [3, 1], [2, 1, 1]), ([3, 1], [2, 1, 1], [3, 1]),
                           (True, False, False), [1, 1, 1, 1, 0], 4)
        self.cgo_book_3 = (([3, 1], [3, 1], [2, 1, 1]), ([2, 2], [2, 2], [2, 2]),
                           (True, False, False), [0, 1, 0, 1, 0], 4)
        self.cgo_book_4 = (([2, 2], [2, 2], [2, 2]), ([2, 2], [2, 2], [2, 2]),
                           (True, False, False), [1, 0, 1, 0, 1], 4)
        self.cgo_book_5 = (([4, 1], [4, 1], [2, 1, 1, 1]), ([4, 1], [2, 1, 1, 1], [4, 1]),
                           (True, False, False), [1, 1, 1, 1, 0, 0, 0], 5)
        self.cgo_book_6 = (([4, 1], [4, 1], [2, 1, 1, 1]), ([3, 2], [2, 2, 1], [3, 2]),
                           (True, False, False), [0, 1, 1, 1, 1, 0, 0], 5)
        self.cgo_book_7 = (([4, 1], [4, 1], [2, 1, 1, 1]), ([3, 1, 1], [3, 1, 1], [3, 1, 1]),
                           (True, False, False), [0, 1, 1, 1, 1, 1, 0], 5)
        self.cgo_book_8 = (([3, 2], [3, 2], [2, 2, 1]), ([3, 2], [2, 2, 1], [3, 2]),
                           (True, False, False), [1, 1, 1, 1, 1, 1, 0], 5)
        self.cgo_book_9 = (([3, 2], [3, 2], [2, 2, 1]), ([3, 1, 1], [3, 1, 1], [3, 1, 1]),
                           (True, False, False), [0, 1, 1, 2, 1, 1, 0], 5)
        self.cgo_book_10 = (([3, 1, 1], [3, 1, 1], [3, 1, 1]), ([3, 1, 1], [3, 1, 1], [3, 1, 1]),
                            (True, False, False), [1, 1, 2, 1, 2, 1, 1], 5)
        self.cgo_book_num_list = list(range(1, 10 + 1))  # [1, ..., 10]

        _, self.s_n_finish_file_name = get_cg_series_finish_s_n_name()
        _, self.s_n_finish_full_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
        self.create_time_dict = {}  # 用于检查计算好的部分不会重复计算

    def teardown_class(self):
        self.protector.protector_teardown()
        pass

    # start with 0xx tests need test by order

    def test_001_create_cg_series(self):
        """for s_n=1, there is no finish db"""
        # check with no db
        for ex in [".pkl", ".txt"]:
            assert not os.path.exists(self.s_n_finish_full_file_name + ex)
        flag, cg_series = load_cg_series(1, [1], [1], is_flag_true_if_not_s_n=True)
        assert flag
        assert cg_series is False
        flag, cg_series = load_cg_series(1, [1], [1], is_flag_true_if_not_s_n=False)
        assert not flag
        assert isinstance(cg_series, str)
        flag, finish_s_n = get_cg_series_finish_s_n()
        assert flag
        assert finish_s_n == 0

        # check create_cg_series_s_n_1
        flag, finish_s_n = create_cg_series(1)
        assert flag
        assert finish_s_n == 1

        # check answer
        for i in self.cgo_case_num_list:
            s_n, yd_1, yd_2, cgo_list = eval("self.cgo_case_{}".format(i))
            if s_n > finish_s_n:
                continue
            _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
            flag, data = CGSeriesInfo(1).query_by_file_name(file_name)
            assert flag
            assert isinstance(data.get("create_time"), str)
            self.create_time_dict["S1"] = data.get("create_time")
            _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
            _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
            assert flag
            assert all(cg_series == np.array(cgo_list))

        for i in self.cgo_book_num_list:
            yd_1_tuple, yd_2_tuple, bool_tuple, cgo_list_tmp, s_n = eval("self.cgo_book_{}".format(i))
            if s_n > finish_s_n:
                continue
            for yd_1, yd_2, bool_value in zip(yd_1_tuple, yd_2_tuple, bool_tuple):
                if bool_value:
                    cgo_list = cgo_list_tmp
                else:
                    cgo_list = reversed(cgo_list_tmp)
                _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
                flag, data = CGSeriesInfo(1).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
                _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
                for ex in [".pkl", ".txt"]:
                    assert os.path.exists(full_file_name + ex)
                    assert os.path.exists(full_finish_file_name + ex)
                flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
                assert flag
                assert all(cg_series == np.array(cgo_list))

        # check finish s_n
        flag, finish_s_n = get_cg_series_finish_s_n()
        assert flag
        assert finish_s_n == 1

        # history_times
        _, finish_file_name = get_cg_series_finish_s_n_name()
        flag, data = CGSeriesInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == np.array([0])
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        assert isinstance(data.get("flags").get("history_times").get("S1"), int)
        assert 0 <= data.get("flags").get("history_times").get("S1") <= 1
        assert data.get("flags").get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"
        flag, data_txt = CGSeriesInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        assert isinstance(data.get("history_times").get("S{}".format(1)), int)
        assert data.get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"

    def test_002_create_cg_series_s_n_2_to_4(self):
        # check create_cg_series_s_n_2_to_4
        flag, finish_s_n = create_cg_series(4)
        assert flag
        assert finish_s_n == 4

        _, file_name = get_cg_series_file_name(1, [1], [1])
        flag, data = CGSeriesInfo(1).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get("S1") == data.get("create_time")

        # check answer
        for i in self.cgo_case_num_list:
            s_n, yd_1, yd_2, cgo_list = eval("self.cgo_case_{}".format(i))
            if s_n > finish_s_n:
                continue
            _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
            flag, data = CGSeriesInfo(s_n).query_by_file_name(file_name)
            assert flag
            # logger.warning("@@@@ {}".format(eval("self.cgo_case_{}".format(i))))
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
            _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                if ex == ".txt" and s_n > 7:
                    continue
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
            assert flag
            assert all(cg_series == np.array(cgo_list))

        for i in self.cgo_book_num_list:
            yd_1_tuple, yd_2_tuple, bool_tuple, cgo_list_tmp, s_n = eval("self.cgo_book_{}".format(i))
            if s_n > finish_s_n:
                continue
            for yd_1, yd_2, bool_value in zip(yd_1_tuple, yd_2_tuple, bool_tuple):
                if bool_value:
                    cgo_list = cgo_list_tmp
                else:
                    cgo_list = [i for i in reversed(cgo_list_tmp)]  # TODO 反向迭代器造出来的和普通生成的list略有不同
                # logger.warning("@@@@ s_n={}, yd_1={}, yd_2={}, bool_value={}, cgo_list_tmp={}, cgo_list={}".format(
                #     s_n, yd_1, yd_2, bool_value, cgo_list_tmp, cgo_list))
                _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
                flag, data = CGSeriesInfo(s_n).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
                _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
                for ex in [".pkl", ".txt"]:
                    if ex == ".txt" and s_n > 7:
                        continue
                    assert os.path.exists(full_file_name + ex)
                    assert os.path.exists(full_finish_file_name + ex)
                flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
                assert flag
                assert all(cg_series == np.array(cgo_list)), "cg_series={}, answer_list={}".format(cg_series, cgo_list)

        # check finish s_n
        flag, finish_s_n = get_cg_series_finish_s_n()
        assert flag
        assert finish_s_n == 4

        # history_times
        _, finish_file_name = get_cg_series_finish_s_n_name()
        flag, data = CGSeriesInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == np.array([0])
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(1, 4 + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        assert data.get("flags").get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"
        flag, data_txt = CGSeriesInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(1, 4 + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        assert data.get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"

    def test_003_create_cg_series_s_n_5_to_9(self):
        # 准备前文
        flag, msg = create_young_diagrams(9)
        assert flag
        assert msg == 9

        flag, msg = create_characters_and_gi(9)
        assert flag
        assert msg == 9

        # check create_cg_series_s_n_5_to_9
        flag, finish_s_n = create_cg_series(9)
        assert flag
        assert finish_s_n == 9

        _, file_name = get_cg_series_file_name(1, [1], [1])
        flag, data = CGSeriesInfo(1).query_by_file_name(file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert self.create_time_dict.get("S1") == data.get("create_time")

        # check answer
        for i in self.cgo_case_num_list:
            s_n, yd_1, yd_2, cgo_list = eval("self.cgo_case_{}".format(i))
            if s_n > finish_s_n:
                continue
            _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
            flag, data = CGSeriesInfo(s_n).query_by_file_name(file_name)
            assert flag
            # logger.warning("@@@@ {}".format(eval("self.cgo_case_{}".format(i))))
            assert isinstance(data.get("create_time"), str)
            _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
            _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
            for ex in [".pkl", ".txt"]:
                if ex == ".txt" and s_n > 7:
                    continue
                assert os.path.exists(full_file_name + ex)
                assert os.path.exists(full_finish_file_name + ex)
            flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
            assert flag
            assert all(cg_series == np.array(cgo_list))

        for i in self.cgo_book_num_list:
            yd_1_tuple, yd_2_tuple, bool_tuple, cgo_list_tmp, s_n = eval("self.cgo_book_{}".format(i))
            if s_n > finish_s_n:
                continue
            for yd_1, yd_2, bool_value in zip(yd_1_tuple, yd_2_tuple, bool_tuple):
                if bool_value:
                    cgo_list = cgo_list_tmp
                else:
                    cgo_list = [i for i in reversed(cgo_list_tmp)]  # TODO 反向迭代器造出来的和普通生成的list略有不同
                _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
                flag, data = CGSeriesInfo(s_n).query_by_file_name(file_name)
                assert flag
                assert isinstance(data.get("create_time"), str)
                _, full_file_name = get_cg_series_file_name(s_n, yd_1, yd_2, is_full_path=True)
                _, full_finish_file_name = get_cg_series_finish_s_n_name(is_full_path=True)
                for ex in [".pkl", ".txt"]:
                    if ex == ".txt" and s_n > 7:
                        continue
                    assert os.path.exists(full_file_name + ex)
                    assert os.path.exists(full_finish_file_name + ex)
                flag, cg_series = load_cg_series(s_n, yd_1, yd_2, is_flag_true_if_not_s_n=True)
                assert flag
                assert all(cg_series == np.array(cgo_list))

        # check finish s_n
        flag, finish_s_n = get_cg_series_finish_s_n()
        assert flag
        assert finish_s_n == 9

        # history_times
        _, finish_file_name = get_cg_series_finish_s_n_name()
        flag, data = CGSeriesInfo(0).query_by_file_name(finish_file_name)
        assert flag
        assert data.get("data") == np.array([0])
        assert isinstance(data.get("flags"), dict)
        assert isinstance(data.get("flags").get("history_times"), dict)
        for i in range(1, 9 + 1):
            assert isinstance(data.get("flags").get("history_times").get("S{}".format(i)), int)
            assert 0 <= data.get("flags").get("history_times").get("S{}".format(i)) <= 1
        assert data.get("flags").get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"
        flag, data_txt = CGSeriesInfo(0).query_txt_by_file_name(finish_file_name)
        assert isinstance(data_txt, str)
        data = eval(data_txt)
        assert isinstance(data, dict)
        for i in range(1, 9 + 1):
            assert isinstance(data.get("history_times").get("S{}".format(i)), int)
        assert data.get("young_diagram_index") == "young diagram list of Sn by young-yamanouchi"
