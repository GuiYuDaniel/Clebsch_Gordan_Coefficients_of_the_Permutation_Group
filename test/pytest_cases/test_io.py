# -*- coding:utf8 -*-
"""
测试src/utils/io功能是否正确执行
"""


import os
import pickle
import pytest
from utils.io import Path, Save, Load, Delete


class TestSave(object):
    """
    test class Save functions
    """

    def setup_class(self):
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.file_path = os.path.join(self.fake_path, "test_results", "fake_cgc", "fake_rst.pkl")
        self.test_data = {"peace": "love"}
        if os.path.exists(self.file_path):
            os.remove(self.file_path)

    def teardown_class(self):
        if os.path.exists(self.file_path):
            os.remove(self.file_path)

    def test_pickle_save_load_and_delete(self):
        # save
        Save.save_pickle(self.test_data, self.file_path)
        assert os.path.exists(self.file_path)
        with open(self.file_path, "rb") as f:
            data = pickle.load(f)
        assert data == self.test_data
        # load
        flag, data = Load.load_pickle(self.file_path)
        assert flag
        assert data == self.test_data
        # delete
        flag, msg = Delete.delete_pickle(self.file_path)
        assert flag
        assert not os.path.exists(self.file_path)
