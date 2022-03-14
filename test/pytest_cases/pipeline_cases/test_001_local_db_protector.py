# -*- coding:utf8 -*-
"""
测试src/db/local_db_protector.py下的功能
"""


import os
import shutil
from db.local_db_protector import DBProtector
from utils.config import singleton_config
from utils.io import Path
from utils.log import get_logger


logger = get_logger(__name__)


class TestDBProtectorInFunc(object):

    def setup(self):
        # create test db
        self.top_path = singleton_config.top_path
        self.tmp_1_path = os.path.join(self.top_path, "tmp_1")
        self.tmp_2_path = os.path.join(self.tmp_1_path, "tmp_2")
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_results_folder = os.path.join(self.fake_path, "test_results")
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipetask_info")
        self.ex_name = ".protected"
        self.ppn_path = os.path.join(self.tmp_2_path, "fake_pipenode_info")
        self.ppl_path = os.path.join(self.tmp_2_path, "fake_pipeline_info")
        self.ppt_path = os.path.join(self.tmp_2_path, "fake_pipetask_info")
        # new tmp folder for test
        os.makedirs(self.tmp_2_path)
        shutil.copytree(self.fake_pipenode_pkl_folder, self.ppn_path)
        shutil.copytree(self.fake_pipeline_pkl_folder, self.ppl_path)
        shutil.copytree(self.fake_pipetask_pkl_folder, self.ppt_path)

    def teardown(self):
        # delete test db
        shutil.rmtree(self.tmp_1_path)

    def test_protector_mv(self):
        protector = DBProtector("tmp_1")

        protector.protector_setup()
        assert os.path.exists(self.tmp_1_path)  # mv新建一个同名目录
        assert not os.path.exists(self.tmp_2_path)  # mv不复制
        protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert os.path.exists(protected_tmp_1_path)  # 被mv保护的目录
        protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert not os.path.exists(protected_tmp_2_path)  # 里层目录跟随，而不是递归mv
        test_file_path = os.path.join(self.tmp_1_path, "test_folder")
        os.makedirs(test_file_path)
        assert os.path.exists(test_file_path)  # 保护期间写入一个文件

        protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)
        assert not os.path.exists(protected_tmp_1_path)
        assert not os.path.exists(protected_tmp_2_path)
        assert not os.path.exists(test_file_path)  # 保护期间写入的文件不应该存在于生产目录

    def test_protector_copy(self):
        protector = DBProtector(os.path.join("tmp_1", "tmp_2"), is_copy=True)

        protector.protector_setup()
        assert os.path.exists(self.tmp_1_path)  # 外层目录
        protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert not os.path.exists(protected_tmp_1_path)  # 外层目录不该被保护
        assert os.path.exists(self.tmp_2_path)  # cp的同名目录
        protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert os.path.exists(protected_tmp_2_path)  # 被保护的目录
        assert os.path.exists(self.ppn_path)  # cp的目录里要啥有啥
        assert os.path.exists(self.ppl_path)
        assert os.path.exists(self.ppt_path)
        protected_ppn_path = self.ppn_path.replace("tmp_2", "tmp_2" + self.ex_name)
        protected_ppl_path = self.ppl_path.replace("tmp_2", "tmp_2" + self.ex_name)
        protected_ppt_path = self.ppt_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert os.path.exists(protected_ppn_path)  # 被保护的目录里也要啥有啥
        assert os.path.exists(protected_ppl_path)
        assert os.path.exists(protected_ppt_path)
        test_file_path = os.path.join(self.tmp_2_path, "test_folder")
        os.makedirs(test_file_path)
        assert os.path.exists(test_file_path)  # 保护期间写入一个文件

        protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)
        assert not os.path.exists(protected_tmp_1_path)
        assert not os.path.exists(protected_tmp_2_path)
        assert not os.path.exists(test_file_path)  # 保护期间写入的文件不应该存在于生产目录
        assert os.path.exists(self.ppn_path)  # 原来的东西不能动
        assert os.path.exists(self.ppl_path)
        assert os.path.exists(self.ppt_path)


class TestDBProtectorInClassMV(object):

    def setup(self):
        # create test db
        self.top_path = singleton_config.top_path
        self.tmp_1_path = os.path.join(self.top_path, "tmp_1")
        self.tmp_2_path = os.path.join(self.tmp_1_path, "tmp_2")
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_results_folder = os.path.join(self.fake_path, "test_results")
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipetask_info")
        self.ex_name = ".protected"
        self.ppn_path = os.path.join(self.tmp_2_path, "fake_pipenode_info")
        self.ppl_path = os.path.join(self.tmp_2_path, "fake_pipeline_info")
        self.ppt_path = os.path.join(self.tmp_2_path, "fake_pipetask_info")
        # new tmp folder for test
        os.makedirs(self.tmp_2_path)
        shutil.copytree(self.fake_pipenode_pkl_folder, self.ppn_path)
        shutil.copytree(self.fake_pipeline_pkl_folder, self.ppl_path)
        shutil.copytree(self.fake_pipetask_pkl_folder, self.ppt_path)

        # protect
        self.protector = DBProtector("tmp_1")
        self.protector.protector_setup()

        # test protector setup
        assert os.path.exists(self.tmp_1_path)  # mv新建一个同名目录
        assert not os.path.exists(self.tmp_2_path)  # mv不复制
        self.protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert os.path.exists(self.protected_tmp_1_path)  # 被mv保护的目录
        self.protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert not os.path.exists(self.protected_tmp_2_path)  # 里层目录跟随，而不是递归mv

    def teardown(self):
        # test protector teardown
        self.protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)
        assert not os.path.exists(self.protected_tmp_1_path)
        assert not os.path.exists(self.protected_tmp_2_path)
        assert not os.path.exists(self.test_file_path)  # 保护期间写入的文件不应该存在于生产目录

        # delete test db
        shutil.rmtree(self.tmp_1_path)

    def test_makedir(self):
        self.test_file_path = os.path.join(self.tmp_1_path, "test_folder")
        os.makedirs(self.test_file_path)
        assert os.path.exists(self.test_file_path)  # 保护期间写入一个文件

    # def test_raise(self):
    #     self.test_file_path = os.path.join(self.tmp_1_path, "test_folder")
    #     os.makedirs(self.test_file_path)
    #     assert os.path.exists(self.test_file_path)  # 保护期间写入一个文件
    #     raise Exception("Raise id support here!")


class TestDBProtectorInClassCP(object):

    def setup(self):
        # create test db
        self.top_path = singleton_config.top_path
        self.tmp_1_path = os.path.join(self.top_path, "tmp_1")
        self.tmp_2_path = os.path.join(self.tmp_1_path, "tmp_2")
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_results_folder = os.path.join(self.fake_path, "test_results")
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipetask_info")
        self.ex_name = ".protected"
        self.ppn_path = os.path.join(self.tmp_2_path, "fake_pipenode_info")
        self.ppl_path = os.path.join(self.tmp_2_path, "fake_pipeline_info")
        self.ppt_path = os.path.join(self.tmp_2_path, "fake_pipetask_info")
        # new tmp folder for test
        os.makedirs(self.tmp_2_path)
        shutil.copytree(self.fake_pipenode_pkl_folder, self.ppn_path)
        shutil.copytree(self.fake_pipeline_pkl_folder, self.ppl_path)
        shutil.copytree(self.fake_pipetask_pkl_folder, self.ppt_path)

        # protect
        self.protector = DBProtector(os.path.join("tmp_1", "tmp_2"), is_copy=True)

        # test protector setup
        self.protector.protector_setup()
        assert os.path.exists(self.tmp_1_path)  # 外层目录
        self.protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert not os.path.exists(self.protected_tmp_1_path)  # 外层目录不该被保护
        assert os.path.exists(self.tmp_2_path)  # cp的同名目录
        self.protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert os.path.exists(self.protected_tmp_2_path)  # 被保护的目录
        assert os.path.exists(self.ppn_path)  # cp的目录里要啥有啥
        assert os.path.exists(self.ppl_path)
        assert os.path.exists(self.ppt_path)
        protected_ppn_path = self.ppn_path.replace("tmp_2", "tmp_2" + self.ex_name)
        protected_ppl_path = self.ppl_path.replace("tmp_2", "tmp_2" + self.ex_name)
        protected_ppt_path = self.ppt_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert os.path.exists(protected_ppn_path)  # 被保护的目录里也要啥有啥
        assert os.path.exists(protected_ppl_path)
        assert os.path.exists(protected_ppt_path)

    def teardown(self):
        self.protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)
        assert not os.path.exists(self.protected_tmp_1_path)
        assert not os.path.exists(self.protected_tmp_2_path)
        assert not os.path.exists(self.test_file_path)  # 保护期间写入的文件不应该存在于生产目录
        assert os.path.exists(self.ppn_path)  # 原来的东西不能动
        assert os.path.exists(self.ppl_path)
        assert os.path.exists(self.ppt_path)

        # delete test db
        shutil.rmtree(self.tmp_1_path)

    def test_makedir(self):
        self.test_file_path = os.path.join(self.tmp_2_path, "test_folder")
        os.makedirs(self.test_file_path)
        assert os.path.exists(self.test_file_path)  # 保护期间写入一个文件

    # def test_makedir_wrong(self):
    #     self.test_file_path = os.path.join(self.tmp_1_path, "test_folder")
    #     os.makedirs(self.test_file_path)
    #     assert os.path.exists(self.test_file_path)  # 保护期间写入一个文件

    # def test_raise(self):
    #     self.test_file_path = os.path.join(self.tmp_2_path, "test_folder")
    #     os.makedirs(self.test_file_path)
    #     assert os.path.exists(self.test_file_path)  # 保护期间写入一个文件
    #     raise Exception("Raise id support here!")


class TestDBProtectorInClassMVWithNoDB(object):

    def setup(self):
        # create test db
        self.top_path = singleton_config.top_path
        self.tmp_1_path = os.path.join(self.top_path, "tmp_1")
        self.tmp_2_path = os.path.join(self.tmp_1_path, "tmp_2")
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_results_folder = os.path.join(self.fake_path, "test_results")
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipetask_info")
        self.ex_name = ".protected"
        self.ppn_path = os.path.join(self.tmp_2_path, "fake_pipenode_info")
        self.ppl_path = os.path.join(self.tmp_2_path, "fake_pipeline_info")
        self.ppt_path = os.path.join(self.tmp_2_path, "fake_pipetask_info")
        # new tmp folder for test
        os.makedirs(self.tmp_1_path)

        # protect
        self.protector = DBProtector(os.path.join("tmp_1", "tmp_2"))  # no tmp_2
        self.protector.protector_setup()

        # test protector setup
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)  # 新建的
        self.protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert not os.path.exists(self.protected_tmp_1_path)
        self.protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert not os.path.exists(self.protected_tmp_2_path)  # 没有原db供备份

    def teardown(self):
        # test protector teardown
        self.protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert not os.path.exists(self.tmp_2_path)
        assert not os.path.exists(self.protected_tmp_1_path)
        assert not os.path.exists(self.protected_tmp_2_path)
        assert not os.path.exists(self.test_file_path_2)  # 保护期间写入的文件不应该存在于生产目录
        assert os.path.exists(self.test_file_path_1)  # 写在外部的还是要有

        # delete test db
        shutil.rmtree(self.tmp_1_path)

    def test_makedir(self):
        self.test_file_path_1 = os.path.join(self.tmp_1_path, "test_folder")
        os.makedirs(self.test_file_path_1)
        assert os.path.exists(self.test_file_path_1)  # 写在外部的
        self.test_file_path_2 = os.path.join(self.tmp_2_path, "test_folder")
        os.makedirs(self.test_file_path_2)
        assert os.path.exists(self.test_file_path_2)  # 保护期间写入一个文件

    # def test_raise(self):
    #     self.test_file_path_1 = os.path.join(self.tmp_1_path, "test_folder")
    #     os.makedirs(self.test_file_path_1)
    #     assert os.path.exists(self.test_file_path_1)  # 写在外部的
    #     self.test_file_path_2 = os.path.join(self.tmp_2_path, "test_folder")
    #     os.makedirs(self.test_file_path_2)
    #     assert os.path.exists(self.test_file_path_2)  # 保护期间写入一个文件
    #     raise Exception("Raise id support here!")
