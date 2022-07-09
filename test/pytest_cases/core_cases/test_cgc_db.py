# -*- coding:utf8 -*-
"""
测试
core/cgc_utils/cgc_local_db.py
core/cgc_utils/cgc_db_typing.py
下所有功能是否正确执行
"""


import copy
import os
import pytest
import shutil
import time
import numpy as np
from conf.cgc_config import top_path, cgc_rst_folder, default_s_n
from core.cgc_utils.cgc_local_db import CGCLocalDb
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo, BranchingLawInfo, YoungTableInfo
from core.cgc_utils.cgc_db_typing import YamanouchiMatrixInfo, CharacterAndGiInfo, CGSeriesInfo
from core.cgc_utils.cgc_db_typing import EigenvaluesInfo, ISFInfo, CGCInfo
from db.local_db_protector import DBProtector
from utils.log import get_logger
# from utils.utils import asctime_2_time


logger = get_logger(__name__)


# TODO delete all pytest.mark.skip("pass")
class TmpTyping(CGCLocalDb):

    def __init__(self, s_n):
        super(TmpTyping, self).__init__()
        self.table_type = "tmp_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": None,
        })
        self.s_n = s_n
        self.txt_limit = 2
        self._init_cgc_static_db_folder()


@pytest.mark.skip("pass")
class TestCGCLocalDb(object):
    """
    test class CGCLocalDb functions
    """
    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

    def setup(self):
        self.db_static_folder = os.path.join(top_path, cgc_rst_folder, "tmp_info")
        if os.path.exists(self.db_static_folder):
            shutil.rmtree(self.db_static_folder)
        self.file_name = "S_3/sigma[1, 1, 1]_mu[2, 1]"  # 故意不带.pkl的
        self.file_name_fake_1 = "S_3"
        self.file_name_fake_2 = "sigma[1, 1, 1]_mu[2, 1]"
        self.test_data = {
            "file_name": self.file_name,
            "data": {1: 1, "2": "2"},
            "flags": None
        }
        self.file_name_pkl = os.path.join(self.db_static_folder, self.file_name) + ".pkl"
        self.file_name_txt = os.path.join(self.db_static_folder, self.file_name) + ".txt"
        self.test_data_same = copy.deepcopy(self.test_data)
        self.update_partial_data = {
            "data": {1: 1, "3": "3"},
            "flags": {"computing_time": "2"}
        }
        self.updated_data = {
            "file_name": self.file_name,
            "data": {1: 1, "3": "3"},  # 不是{1: 1, "2": "2", "3": "3"}哦
            "flags": {"computing_time": "2"}
        }  # 这里没写时间项

    def teardown(self):
        if os.path.exists(self.db_static_folder):
            shutil.rmtree(self.db_static_folder)

    def teardown_class(self):
        self.protector.protector_teardown()

    def test_001_init_cgc_static_db_folder(self):  # 001使得执行顺序靠前
        assert not os.path.exists(self.db_static_folder)
        db_info = TmpTyping(s_n=1)
        assert os.path.exists(self.db_static_folder)

    def test_db_op(self):
        db_info = TmpTyping(s_n=-1)
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        # test insert db
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        # test query db
        flag, data = db_info.query_by_file_name(self.file_name)
        assert flag
        assert data is not False
        for key in self.test_data:
            assert key in data
            assert self.test_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        assert os.path.exists(self.file_name_pkl)

        flag, data = db_info.query_by_file_name(self.file_name_fake_1)
        assert flag
        assert data is False
        flag, data = db_info.query_by_file_name(self.file_name_fake_2)
        assert flag
        assert data is False

        # test update
        flag, msg = db_info.update_by_file_name(self.file_name_fake_1, self.update_partial_data)
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False

        m_time_before = os.stat(self.file_name_pkl).st_mtime
        flag, msg = db_info.update_by_file_name(self.file_name, self.test_data)
        assert flag
        assert msg
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before == m_time_after

        m_time_before = os.stat(self.file_name_pkl).st_mtime
        flag, msg = db_info.update_by_file_name(self.file_name, self.update_partial_data)  # 真的
        assert flag
        assert msg
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before != m_time_after
        flag, data = db_info.query_by_file_name(self.file_name)
        assert flag
        assert data is not False
        for key in self.updated_data:
            assert key in data
            assert self.updated_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        assert m_time_after == os.stat(self.file_name_pkl).st_mtime

        m_time_before = os.stat(self.file_name_pkl).st_mtime  # update again
        flag, msg = db_info.update_by_file_name(self.file_name, self.update_partial_data)
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before == m_time_after
        flag, data = db_info.query_by_file_name(self.file_name)
        assert flag
        assert data is not False
        for key in self.updated_data:
            assert key in data
            assert self.updated_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        assert m_time_after == os.stat(self.file_name_pkl).st_mtime
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        # test delete
        flag, msg = db_info.delete_by_file_name(self.file_name_fake_1)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)

        flag, msg = db_info.delete({"file_name": self.file_name_fake_2})  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)

        flag, msg = db_info.delete_by_file_name(self.test_data.get("file_name"))  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_by_file_name(self.test_data.get("file_name"))  # delete again
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

    def test_within_txt_limit(self):
        db_info = TmpTyping(s_n=2)  # s_n = txt_limit, auto判断会存
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        # test insert db 混入项 得保证txt不会动db
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        m_pkl_time_before = os.stat(self.file_name_pkl).st_mtime

        # test insert txt
        flag, msg = db_info.insert_txt(self.test_data)
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_txt)
        m_pkl_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_pkl_time_before == m_pkl_time_after

        # insert txt again
        m_time_before = os.stat(self.file_name_txt).st_mtime
        flag, msg = db_info.insert_txt(self.test_data)
        logger.info("Error is supposed here!")
        assert self.test_data == self.test_data_same
        assert flag is False
        assert msg
        assert os.path.exists(self.file_name_txt)
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before == m_time_after

        # test query txt
        flag, data = db_info.query_txt_by_file_name(self.file_name)
        assert flag
        assert data == str(self.test_data_same.get("data"))

        flag, data = db_info.query_txt_by_file_name(self.file_name_fake_1)
        assert flag
        assert data is False
        flag, data = db_info.query_txt_by_file_name(self.file_name_fake_2)
        assert flag
        assert data is False

        # test update txt
        m_time_before = os.stat(self.file_name_txt).st_mtime
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.test_data)
        assert flag
        assert msg
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before == m_time_after

        m_time_before = os.stat(self.file_name_txt).st_mtime
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.update_partial_data)
        assert flag
        assert msg
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before != m_time_after
        flag, data = db_info.query_txt_by_file_name(self.file_name)
        assert flag
        assert data == str(self.updated_data.get("data"))
        assert m_time_after == os.stat(self.file_name_txt).st_mtime

        m_time_before = os.stat(self.file_name_txt).st_mtime  # update again
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.update_partial_data)
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before == m_time_after
        flag, data = db_info.query_txt_by_file_name(self.file_name)
        assert flag
        assert data == str(self.updated_data.get("data"))
        assert m_time_after == os.stat(self.file_name_txt).st_mtime

        m_pkl_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_pkl_time_before == m_pkl_time_after

        # test delete txt
        flag, msg = db_info.delete_txt_by_file_name(self.file_name_fake_2)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.test_data.get("file_name"))  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.test_data.get("file_name"))  # delete again
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        m_pkl_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_pkl_time_before == m_pkl_time_after

        # delete db
        flag, msg = db_info.delete_by_file_name(self.test_data.get("file_name"))  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

    def test_out_of_txt_limit(self):
        db_info = TmpTyping(s_n=3)  # s_n > txt_limit, auto判断结果为不存
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

        # test insert txt
        flag, msg = db_info.insert_txt(self.test_data)
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is False
        assert not os.path.exists(self.file_name_txt)

        flag, msg = db_info.insert_txt(self.test_data, is_auto_limit=False)
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_txt)

        # test query txt
        flag, data = db_info.query_txt_by_file_name(self.file_name)
        assert flag
        assert data is False

        flag, data = db_info.query_txt_by_file_name(self.file_name, is_auto_limit=False)
        assert flag
        assert data == str(self.test_data_same.get("data"))

        flag, data = db_info.query_txt_by_file_name(self.file_name_fake_1, is_auto_limit=False)
        assert flag
        assert data is False

        # test update txt
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.update_partial_data)
        assert flag
        assert msg is False

        flag, msg = db_info.update_txt_by_file_name(self.file_name_fake_1, self.update_partial_data)
        assert flag
        assert msg is False

        m_time_before = os.stat(self.file_name_txt).st_mtime  # data一致，无需实质更新
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.test_data, is_auto_limit=False)
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before == m_time_after

        m_time_before = os.stat(self.file_name_txt).st_mtime  # data不一致，需实质更新
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.update_partial_data, is_auto_limit=False)
        assert flag
        assert msg
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before != m_time_after
        flag, data = db_info.query_txt_by_file_name(self.file_name, is_auto_limit=False)
        assert flag
        assert data is not False
        assert data == str(self.updated_data.get("data"))
        assert m_time_after == os.stat(self.file_name_txt).st_mtime

        m_time_before = os.stat(self.file_name_txt).st_mtime  # update again
        flag, msg = db_info.update_txt_by_file_name(self.file_name, self.update_partial_data, is_auto_limit=False)
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_txt).st_mtime
        assert m_time_before == m_time_after

        # test delete txt
        flag, msg = db_info.delete_txt_by_file_name(self.file_name_fake_2)  # 应该是找不到的
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.file_name_fake_2, is_auto_limit=False)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.test_data.get("file_name"))  # 不可以
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.test_data.get("file_name"), is_auto_limit=False)  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_txt)

        flag, msg = db_info.delete_txt_by_file_name(self.test_data.get("file_name"), is_auto_limit=False)  # delete again
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)


@pytest.mark.skip("pass")
class TestYoungDiagramInfo(object):
    """
    test python file cgc_db_typing.py:YoungDiagramInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = "S{}".format(self.fake_finish_s_n)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": [[3], [2, 1], [1, 1, 1]],
            "flags": {"speed_time": 0.1}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]],
                                   "flags": {"speed_time": 0.2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": [],
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 0.01}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = {"flags": {"finish_s_n": 1001,
                                                        "history_times": {"S1000": 0.01,
                                                                          "S1001": 0.03}
                                                        }}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_young_diagrams_info(self):
        db_info = YoungDiagramInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_young_diagrams_info_txt(self):
        db_info = YoungDiagramInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_young_diagrams_finish_s_n_info(self):
        db_info = YoungDiagramInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_young_diagrams_finish_s_n_info_txt(self):
        db_info = YoungDiagramInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = YoungDiagramInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = YoungDiagramInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


@pytest.mark.skip("pass")
class TestBranchingLawInfo(object):
    """
    test python file cgc_db_typing.py:BranchingLawInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_nu = [2, 1]
        self.fake_file_name = "S{}/{}".format(self.fake_finish_s_n, self.fake_nu)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": {
                "BL_num": 2,
                "rows": [1, 0],
                "cols": [0, 1],
                "before_YD": [[2], [1, 1]]
            },
            "flags": {"speed_time": 0.1}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": {"fake": "fake"},
                                   "flags": {"speed_time": 0.2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": {},
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 0.01}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = {"flags": {"finish_s_n": 1001,
                                                        "history_times": {"S1000": 0.01,
                                                                          "S1001": 0.03}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_branching_laws_info(self):
        db_info = BranchingLawInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_branching_laws_info_txt(self):
        db_info = BranchingLawInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_branching_laws_finish_s_n_info(self):
        db_info = BranchingLawInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_branching_laws_finish_s_n_info_txt(self):
        db_info = BranchingLawInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = BranchingLawInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = BranchingLawInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


@pytest.mark.skip("pass")
class TestYoungTableInfo(object):
    """
    test python file cgc_db_typing.py:YoungTableInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_nu = [2, 1]
        self.fake_file_name = "S{}/{}".format(self.fake_finish_s_n, self.fake_nu)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": {
                "1": [[1, 2], [3]],
                "2": [[1, 3], [2]]
            },
            "flags": {"speed_time": 0,
                      "total_num": 2}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": {"1": ["fake"]},
                                   "flags": {"speed_time": 2,
                                             "total_num": 1}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": {},
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = {"flags": {"finish_s_n": 1001,
                                                        "history_times": {"S1000": 1,
                                                                          "S1001": 3}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

        self.fake_yt_num_name = "S{}/{}_num".format(self.fake_finish_s_n, self.fake_nu)
        self.fake_yt_num_table = {
            "file_name": self.fake_yt_num_name,
            "data": {"total_num": 2},
            "flags": {}
        }
        self.fake_yt_num_table_copy = copy.deepcopy(self.fake_yt_num_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_young_tableaux_info(self):
        db_info = YoungTableInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_young_tableaux_info_txt(self):
        db_info = YoungTableInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_young_tableaux_finish_s_n_info(self):
        db_info = YoungTableInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_young_tableaux_finish_s_n_info_txt(self):
        db_info = YoungTableInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = YoungTableInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = YoungTableInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False

    def test_fake_yt_num(self):
        db_info = YoungTableInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_yt_num_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_yt_num_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_yt_num_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_yt_num_name, self.fake_yt_num_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_yt_num_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") == last_write_time_1
        assert data_2.get("data") == self.fake_yt_num_table_copy.get("data")
        assert data_2.get("flags") == self.fake_yt_num_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_yt_num_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_yt_num_name})
        assert flag
        assert data_3 is False


@pytest.mark.skip("pass")
class TestYamanouchiMatrixInfo(object):
    """
    test python file cgc_db_typing.py:YamanouchiMatrixInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_nu = [2, 1]
        self.fake_ij = (2, 3,)
        self.fake_file_name = "S{}/{}/ij{}".format(self.fake_finish_s_n, self.fake_nu, self.fake_ij)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": np.array([[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]]),
            "flags": {"speed_time": 0,
                      "total_num": 2}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": ["fake"],
                                   "flags": {"speed_time": 2,
                                             "total_num": 1}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": [],
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = {"flags": {"finish_s_n": 1001,
                                                        "history_times": {"S1000": 1,
                                                                          "S1001": 3}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_yamanouchi_matrix_info(self):
        db_info = YamanouchiMatrixInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_yamanouchi_matrix_info_txt(self):
        db_info = YamanouchiMatrixInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_yamanouchi_matrix_finish_s_n_info(self):
        db_info = YamanouchiMatrixInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_yamanouchi_matrix_finish_s_n_info_txt(self):
        db_info = YamanouchiMatrixInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = YamanouchiMatrixInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = YamanouchiMatrixInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


@pytest.mark.skip("pass")
class TestCharacterAndGiInfo(object):
    """
    test python file cgc_db_typing.py:CharacterAndGiInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = "S{}".format(self.fake_finish_s_n)
        self.character_2_matrix = np.array([[1, 1], [1, -1]])
        self.gi_2 = np.array([1, 1])
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": {"character": self.character_2_matrix, "gi": self.gi_2},
            "flags": {"speed_time": 0,
                      "young_diagram_index": [[2], [1, 1]]}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": {"fake": "fake"},
                                   "flags": {"speed_time": 2,
                                             "young_diagram_index": [[3], [2, 1], [1, 1, 1]]}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": {},
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = {"flags": {"finish_s_n": 1001,
                                                        "history_times": {"S1000": 1,
                                                                          "S1001": 3}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_character_info(self):
        db_info = CharacterAndGiInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert (data.get("data").get("character") == self.fake_table_copy["data"]["character"]).all()
        assert (data.get("data").get("gi") == self.fake_table_copy["data"]["gi"]).all()
        # assert data == self.fake_table_copy
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_character_info_txt(self):
        db_info = CharacterAndGiInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_character_finish_s_n_info(self):
        db_info = CharacterAndGiInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_character_finish_s_n_info_txt(self):
        db_info = CharacterAndGiInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = CharacterAndGiInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = CharacterAndGiInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


@pytest.mark.skip("pass")
class TestCGSeriesInfo(object):
    """
    test python file cgc_db_typing.py:CGSeriesInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = "S{}/{}_{}".format(self.fake_finish_s_n, [3], [2, 1])
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": np.array([0, 1, 0]),
            "flags": {"speed_time": 0}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": np.array([0, 1, 0, 10000]),
                                   "flags": {"speed_time": 2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": np.array([0]),
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1},
                      "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = \
            {"flags": {"finish_s_n": 1001,
                       "history_times": {"S1000": 1, "S1001": 3},
                       "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_cg_series_info(self):
        db_info = CGSeriesInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert (data.get("data") == self.fake_table_copy["data"]).all()
        assert data.get("flags") == self.fake_table_copy["flags"]
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert (data_2.get("data") == self.fake_partial_table_copy.get("data")).all()
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_cg_series_info_txt(self):
        db_info = CGSeriesInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_cg_series_finish_s_n_info(self):
        db_info = CGSeriesInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_cg_series_finish_s_n_info_txt(self):
        db_info = CGSeriesInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = CGSeriesInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = CGSeriesInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


@pytest.mark.skip("pass")
class TestEigenvaluesInfo(object):
    """
    test python file cgc_db_typing.py:EigenvaluesInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = "S{}".format(self.fake_finish_s_n)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": [15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15],
            "flags": {"speed_time": 0}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": [15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15, 10000],
                                   "flags": {"speed_time": 2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": [],
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1},
                      "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = \
            {"flags": {"finish_s_n": 1001,
                       "history_times": {"S1000": 1, "S1001": 3},
                       "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_eigenvalue_info(self):
        db_info = EigenvaluesInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data.get("data") == self.fake_table_copy["data"]
        assert data.get("flags") == self.fake_table_copy["flags"]
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_eigenvalue_info_txt(self):
        db_info = EigenvaluesInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_eigenvalue_finish_s_n_info(self):
        db_info = EigenvaluesInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_eigenvalue_finish_s_n_info_txt(self):
        db_info = EigenvaluesInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = EigenvaluesInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = EigenvaluesInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


class TestISFInfo(object):
    """
    test python file cgc_db_typing.py:ISFInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = os.path.join("S{}", "{}_{}", "{}\'").format(self.fake_finish_s_n, [2], [2], [1])
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": {"rows": [([1], [1])],
                     "cols": [[2]],
                     "isf": np.array([[1]])},
            "flags": {"speed_time": 0}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": {"rows": [([1], [1])],
                                   "cols": [[1, 1]],
                                   "isf": np.array([[1]])},
                                   "flags": {"speed_time": 2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": {},
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = \
            {"flags": {"finish_s_n": 1001,
                       "history_times": {"S1000": 1, "S1001": 3}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_isf_info(self):
        db_info = ISFInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data.get("data") == self.fake_table_copy["data"]
        assert data.get("flags") == self.fake_table_copy["flags"]
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_isf_info_txt(self):
        db_info = ISFInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_isf_finish_s_n_info(self):
        db_info = ISFInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_isf_finish_s_n_info_txt(self):
        db_info = ISFInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = ISFInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = ISFInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False


class TestCGCInfo(object):
    """
    test python file cgc_db_typing.py:CGCInfo class
    TestCGCLocalDb中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        self.fake_finish_s_n = 1000
        self.fake_file_name = os.path.join("S{}", "{}_{}", "{}_{}_m{}").format(self.fake_finish_s_n,
                                                                               [1], [1], [1], "", 1)
        self.fake_table = {
            "file_name": self.fake_file_name,
            "data": {(1, 1,): 1, "N": 1},
            "flags": {"speed_time": 0}
        }
        self.fake_table_copy = copy.deepcopy(self.fake_table)
        self.fake_partial_table = {"data": {(1, 1,): -100, "N": 100},
                                   "flags": {"speed_time": 2}}
        self.fake_partial_table_copy = copy.deepcopy(self.fake_partial_table)

        self.fake_finish_s_n_name = "Finish_Sn"
        self.fake_finish_s_n_table = {
            "file_name": self.fake_finish_s_n_name,
            "data": {},
            "flags": {"finish_s_n": self.fake_finish_s_n,
                      "history_times": {"S1000": 1}}
        }
        self.fake_finish_s_n_table_copy = copy.deepcopy(self.fake_finish_s_n_table)
        self.fake_finish_s_n_partial_table = \
            {"flags": {"finish_s_n": 1001,
                       "history_times": {"S1000": 1, "S1001": 3}}}
        self.fake_finish_s_n_partial_table_copy = copy.deepcopy(self.fake_finish_s_n_partial_table)

    def teardown(self):
        self.protector.protector_teardown()

    def test_cgc_info(self):
        db_info = CGCInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_file_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data.get("data") == self.fake_table_copy["data"]
        assert data.get("flags") == self.fake_table_copy["flags"]
        # update
        time.sleep(1.1)
        flag, msg = db_info.update_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_file_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("data") == self.fake_partial_table_copy.get("data")
        assert data_2.get("flags") == self.fake_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_cgc_info_txt(self):
        db_info = CGCInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_file_name)
        assert flag
        assert data == str(self.fake_table_copy.get("data"))
        # update
        flag, msg = db_info.update_txt_by_file_name(self.fake_file_name, self.fake_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_file_name})
        assert data_2 == str(self.fake_partial_table_copy.get("data"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_file_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_file_name})
        assert flag
        assert data_3 is False

    def test_cgc_finish_s_n_info(self):
        db_info = CGCInfo(self.fake_finish_s_n)
        # insert
        flag, msg = db_info.insert(self.fake_finish_s_n_table)
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.fake_finish_s_n_table_copy
        # update
        flag, msg = db_info.update_by_file_name(self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("data") == self.fake_finish_s_n_table_copy.get("data")  # not update
        assert data_2.get("flags") == self.fake_finish_s_n_partial_table_copy.get("flags")
        # test delete
        flag, msg = db_info.delete_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info.query({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_cgc_finish_s_n_info_txt(self):
        db_info = CGCInfo(0)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is True
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data == str(self.fake_finish_s_n_table_copy.get("flags"))
        # update
        flag, msg = db_info.update_txt_by_file_name(
            self.fake_finish_s_n_name, self.fake_finish_s_n_partial_table, point_key="flags")
        assert flag
        assert msg is True
        flag, data_2 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert data_2 == str(self.fake_finish_s_n_partial_table_copy.get("flags"))
        # test delete
        flag, msg = db_info.delete_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert msg is True
        flag, data_3 = db_info._query_txt({"file_name": self.fake_finish_s_n_name})
        assert flag
        assert data_3 is False

    def test_limit_txt(self):
        db_info = CGCInfo(default_s_n + 1)
        # insert
        flag, msg = db_info.insert_txt(self.fake_finish_s_n_table, point_key="flags")
        assert flag
        assert msg is False
        # query
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
        # query again
        db_info = CGCInfo(0)
        flag, data = db_info.query_txt_by_file_name(self.fake_finish_s_n_name)
        assert flag
        assert data is False
