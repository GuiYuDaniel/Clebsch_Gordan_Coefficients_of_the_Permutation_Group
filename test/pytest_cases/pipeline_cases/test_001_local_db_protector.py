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
        self.top_path = singleton_config.top_path
        self.tmp_1_path = os.path.join(self.top_path, "tmp_1")
        self.tmp_2_path = os.path.join(self.tmp_1_path, "tmp_2")
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_results_folder = os.path.join(self.fake_path, "test_results")
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_results_folder, "fake_pipetask_info")
        self.ex_name = ".protected"
        # new tmp folder for test
        os.makedirs(self.tmp_2_path)
        shutil.copytree(self.fake_pipenode_pkl_folder, self.tmp_2_path)
        shutil.copytree(self.fake_pipeline_pkl_folder, self.tmp_2_path)
        shutil.copytree(self.fake_pipetask_pkl_folder, self.tmp_2_path)

    def teardown(self):
        shutil.rmtree(self.tmp_1_path)

    def test_protector_mv_1(self):
        protector = DBProtector("tmp_1")

        protector.protector_setup()
        assert os.path.exists(self.tmp_1_path)
        assert not os.path.exists(self.tmp_2_path)
        protected_tmp_1_path = self.tmp_1_path.replace("tmp_1", "tmp_1" + self.ex_name)
        assert os.path.exists(protected_tmp_1_path)
        protected_tmp_2_path = self.tmp_2_path.replace("tmp_2", "tmp_2" + self.ex_name)
        assert os.path.exists(protected_tmp_2_path)
        test_file_path = os.path.join(self.tmp_1_path, "test_folder")
        os.makedirs(test_file_path)
        assert test_file_path

        protector.protector_teardown()
        assert os.path.exists(self.tmp_1_path)
        assert os.path.exists(self.tmp_2_path)
        assert not os.path.exists(protected_tmp_1_path)
        assert not os.path.exists(protected_tmp_2_path)
        assert test_file_path

    def test_protector_copy(self):
        pass


class TestDBProtectorInClass(object):
    pass
