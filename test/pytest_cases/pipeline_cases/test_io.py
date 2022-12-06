# -*- coding:utf8 -*-
"""
测试src/utils/io功能是否正确执行
"""


import os
import pickle
import pytest
from db.local_db_protector import DBProtector
from utils.config import singleton_config
from utils.io import Path, Save, Load, Delete
from utils.io import _make_file_path_with_lawful_file_type


# class TestPath(object):
#     """test class Path functions"""
#
#     # TODO 这个类本质上是使用了os.path.abspath方法，暂时没有找到法二，故暂不检查
#     pass


class TestIOHelper(object):
    """test functions in io.py but not in class"""

    def setup_class(self):
        self.top_path = singleton_config.top_path
        self.root_path = "/"
        self.readme_path = os.path.join(self.top_path, "README.md")
        self.gitignore_path = os.path.join(self.top_path, ".gitignore")

    def teardown_class(self):
        pass

    def test_common_situation(self):
        # file_path没有带扩展名，加上file_type返回
        flag, rst = _make_file_path_with_lawful_file_type(self.top_path, ".pkl")
        assert rst == self.top_path + ".pkl"

        # 带了正确的扩展名，原样返回
        top_path_pkl = self.top_path + ".pkl"
        flag, rst = _make_file_path_with_lawful_file_type(top_path_pkl, ".pkl")
        assert rst == top_path_pkl

        # 带了正确的扩展名，原样返回
        flag, rst = _make_file_path_with_lawful_file_type(self.readme_path, ".md")
        assert rst == self.readme_path

        # file_path没有带扩展名，加上file_type返回
        flag, rst = _make_file_path_with_lawful_file_type(self.gitignore_path, ".gitignore")
        assert rst == self.gitignore_path + ".gitignore"
        assert rst != self.gitignore_path

        # file_path没有带扩展名，加上file_type返回
        flag, rst = _make_file_path_with_lawful_file_type(self.gitignore_path, ".pkl")
        assert rst == self.gitignore_path + ".pkl"
        assert rst != self.gitignore_path

        # 带错扩展名了，false
        flag, rst = _make_file_path_with_lawful_file_type(self.readme_path, ".pkl")
        assert not flag

    def test_rare_situation(self):
        # file_path没有带扩展名，加上file_type返回
        # 即使file_path是一个目录
        flag, rst = _make_file_path_with_lawful_file_type("/guiyu/guiyu/", ".pkl")
        assert rst == "/guiyu/guiyu/.pkl"

        # file_path没有带扩展名，加上file_type返回
        # 即使file_path就是根目录
        flag, rst = _make_file_path_with_lawful_file_type("/", ".pkl")
        assert rst == "/.pkl"

        # file_path没有带扩展名，加上file_type返回
        flag, rst = _make_file_path_with_lawful_file_type("/.", ".txt")
        assert rst == "/..txt"

        # file_path没有带扩展名，加上file_type返回
        flag, rst = _make_file_path_with_lawful_file_type("/..", ".txt")
        assert rst == "/...txt"

        # 带了正确的扩展名，原样返回
        flag, rst = _make_file_path_with_lawful_file_type("/guiyu/guiyu.txt.txt.txt", ".txt")
        assert rst == "/guiyu/guiyu.txt.txt.txt"

    def test_common_error(self):
        # file_type 错误
        flag, rst = _make_file_path_with_lawful_file_type("/guiyu", "pkl")
        assert not flag

        # file_type 错误
        flag, rst = _make_file_path_with_lawful_file_type("/guiyu", 1)
        assert not flag

        # file_type 错误
        flag, rst = _make_file_path_with_lawful_file_type("/guiyu", [])
        assert not flag

        # file_path 错误
        flag, rst = _make_file_path_with_lawful_file_type(1, ".pkl")
        assert not flag

        # file_path 错误
        flag, rst = _make_file_path_with_lawful_file_type([], ".pkl")
        assert not flag

        # file_path 错误
        flag, rst = _make_file_path_with_lawful_file_type({}, ".pkl")
        assert not flag

        # file_path 错误
        flag, rst = _make_file_path_with_lawful_file_type("etc", ".pkl")
        assert not flag

        # file_path 错误
        flag, rst = _make_file_path_with_lawful_file_type("results/guiyu", ".pkl")
        assert not flag

    def test_rare_error(self):
        pass


class TestIO(object):
    """
    test class Save Load Delete functions
    """

    def setup_class(self):
        self.result_folder = singleton_config.result_folder
        self.sys_db_name = singleton_config.sys_db_name

        self.protector = DBProtector(self.sys_db_name, extension_name=".test_db_protected")
        self.protector.protector_setup()

        self.fake_path = os.path.join(self.result_folder, self.sys_db_name, "fake", "test")
        self.file_path_without_ext = os.path.join(self.fake_path, "test_results", "fake_cgc", "fake_rst")
        # 注意，使用save/load/delete _ pkl/txt 这些方法时可以写明扩展名，也可以不写
        self.file_path_pkl = self.file_path_without_ext + ".pkl"
        self.file_path_txt = self.file_path_without_ext + ".txt"
        self.test_data = {"peace": "love"}

    def teardown_class(self):
        self.protector.protector_teardown()

    def test_pickle_save_load_and_delete_pkl(self):
        # save
        flag, msg = Save.save_pickle(self.test_data, self.file_path_without_ext)
        assert flag
        assert msg
        assert os.path.exists(self.file_path_pkl)
        assert not os.path.exists(self.file_path_without_ext)
        assert not os.path.exists(self.file_path_txt)
        with open(self.file_path_pkl, "rb") as f:
            data = pickle.load(f)
        assert data == self.test_data
        # load
        flag, data = Load.load_pickle(self.file_path_pkl)
        assert flag
        assert data == self.test_data
        # delete
        flag, msg = Delete.delete_pickle(self.file_path_without_ext)
        assert flag
        assert not os.path.exists(self.file_path_pkl)

    def test_pickle_save_load_and_delete_txt(self):
        # save
        flag, msg = Save.save_txt(str(self.test_data), self.file_path_txt)
        assert flag
        assert msg
        assert os.path.exists(self.file_path_txt)
        assert not os.path.exists(self.file_path_without_ext)
        assert not os.path.exists(self.file_path_pkl)
        with open(self.file_path_txt, "r") as f:
            data = f.read()
        assert data == str(self.test_data)
        # load
        flag, data = Load.load_txt(self.file_path_without_ext)
        assert flag
        assert data == str(self.test_data)
        # delete
        flag, msg = Delete.delete_txt(self.file_path_txt)
        assert flag
        assert not os.path.exists(self.file_path_txt)
