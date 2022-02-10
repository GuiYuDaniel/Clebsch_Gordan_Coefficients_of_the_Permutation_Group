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
from conf.cgc_config import top_path, cgc_rst_folder
from core.cgc_utils.cgc_local_db import CGCLocalDb
from utils.log import get_logger


logger = get_logger(__name__)


class TmpTyping(CGCLocalDb):

    def __init__(self, s_n):
        super(TmpTyping, self).__init__()
        self.table_type = "tmp_info"
        self.map_id = "file_name"
        self.s_n = s_n
        self.txt_limit = 2
        self._init_cgc_static_db_folder()


class TestCGCLocalDb(object):
    """
    test class CGCLocalDb functions
    """

    def setup_class(self):
        self.db_static_folder = os.path.join(top_path, cgc_rst_folder, "tmp_info")
        if os.path.exists(self.db_static_folder):
            shutil.rmtree(self.db_static_folder)
        self.file_name = "S_3/sigema[1, 1, 1]_miu[2, 1]"
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

    def teardown_class(self):
        if os.path.exists(self.db_static_folder):
            shutil.rmtree(self.db_static_folder)

    def test_001_init_cgc_static_db_folder(self):  # 001使得执行顺序靠前
        assert not os.path.exists(self.db_static_folder)
        db_info = TmpTyping(s_n=2)
        assert os.path.exists(self.db_static_folder)

    def test_insert_update_query_and_delete_with_pkl_and_txt(self):
        db_info = TmpTyping(s_n=2)  # s_n <= txt_limit才会存txt 2 <= 2
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        # test insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)
        assert os.path.exists(self.file_name_txt)
        # test query
        flag, data = db_info.query_by_file_name(self.file_name)
        assert flag
        assert data is not False
        for key in self.test_data:
            assert key in data
            assert self.test_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        flag, data_same = db_info.query({"file_name": self.file_name})
        assert flag
        assert data_same is not False
        assert data == data_same
        assert (data is not data_same)  # 保证它们俩的地址不是一个
        assert os.path.exists(self.file_name_pkl)
        assert os.path.exists(self.file_name_txt)
        flag, data = db_info.query_by_file_name(self.file_name_fake_1)
        assert flag
        assert not data
        # test delete
        flag, msg = db_info.delete_by_file_name(self.file_name_fake_1)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)
        assert os.path.exists(self.file_name_txt)
        flag, msg = db_info.delete({"file_name": self.file_name_fake_2})  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)
        assert os.path.exists(self.file_name_txt)
        flag, msg = db_info.delete_by_file_name(self.test_data.get("file_name"))  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)

    def test_insert_update_query_and_delete_with_pkl_only(self):
        db_info = TmpTyping(s_n=3)  # s_n > txt_limit不会存txt 3 > 2
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        # test insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        # test query
        flag, data = db_info.query_by_file_name(self.file_name)
        assert flag
        assert data is not False
        for key in self.test_data:
            assert key in data
            assert self.test_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        assert not os.path.exists(self.file_name_txt)
        # test delete
        flag, msg = db_info.delete_by_file_name(self.file_name_fake_1)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        flag, msg = db_info.delete({"file_name": self.file_name_fake_2})  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)
        flag, msg = db_info.delete_by_file_name(self.test_data.get("file_name"))  # 可以找到
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)
        assert not os.path.exists(self.file_name_txt)


