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
from conf.cgc_config import top_path, cgc_rst_folder
from core.cgc_utils.cgc_local_db import CGCLocalDb
from utils.log import get_logger
from utils.utils import asctime_2_time


logger = get_logger(__name__)


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


class TestCGCLocalDb(object):
    """
    test class CGCLocalDb functions
    """

    def setup(self):
        self.db_static_folder = os.path.join(top_path, cgc_rst_folder, "tmp_info")
        if os.path.exists(self.db_static_folder):
            shutil.rmtree(self.db_static_folder)
        self.file_name = "S_3/sigema[1, 1, 1]_miu[2, 1]"  # 故意不带.pkl的
        self.file_name_fake_1 = "S_3"
        self.file_name_fake_2 = "sigema[1, 1, 1]_miu[2, 1]"
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
