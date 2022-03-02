# -*- coding:utf8 -*-
"""
测试pipenode能否被成功实例化，检查，调用
"""


import os
import pytest
import shutil
from db.typing import PipeNodeInfo
from utils.log import get_logger
from utils.io import Path
from utils.json_op import Json
from pipeline.pipenode import PipeNode


logger = get_logger(__name__)


class TestPipeNode(object):

    def setup_class(self):
        # 这段是为了准备新建ppn
        self.top_path = Path._get_full_path()
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.whereami_file_path = os.path.join(self.fake_path, "fake_workflow_whereami.json")
        self.workflow_conf_dict = Json.file_to_json(self.whereami_file_path)[0]
        self.workflow_conf_dict["prep_nodes"] = []
        self.answer_ppn_dict_partial = {"pipenode_name": "first_node",
                                        "func_des": ["None.test.fake.fake_core", "where_am_i", ""],
                                        "func_str": "from test.fake.fake_core import where_am_i",
                                        "type": "cold",
                                        "inputs": ["point_path"],
                                        "outputs": ["results"],
                                        "next_nodes": [],
                                        "prep_nodes": [],
                                        "outputs_r": {},
                                        "flags": []}
        # 这段是为了准备load ppn
        self.pipenode_db_folder = os.path.join(self.top_path, "results", "pipenode_info")
        self.pipenode_id = "test_load_pipenode_id"
        self.fake_pipenode_pkl_path = os.path.join(self.fake_path,
                                                   "test_results", "fake_pipenode_info", "test_load_pipenode_id.pkl")
        if os.path.exists(self.pipenode_db_folder):
            all_file_name_list = os.listdir(self.pipenode_db_folder)
            all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
            for test_file in all_test_ahead_list:
                os.remove(os.path.join(self.pipenode_db_folder, test_file))
        else:
            os.makedirs(self.pipenode_db_folder)
        shutil.copy(self.fake_pipenode_pkl_path, self.pipenode_db_folder)
        self.appeared_id_list = []

    def teardown_class(self):
        # 删除粘贴过去的ppn pkl
        all_file_name_list = os.listdir(self.pipenode_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipenode_db_folder, test_file))
        # 防止delete失败，尝试删除所有测试中产生的id
        for appeared_id in self.appeared_id_list:
            file_path = os.path.join(self.pipenode_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                logger.warning("find ppn_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)

    def test_create_blank_and_load_ppn(self):
        self.appeared_id_list.append(self.pipenode_id)
        ppn = PipeNode().load_by_id(self.pipenode_id)
        ppn_none_dict = {k: ppn.__dict__[k] for k in ppn.__dict__ if ppn.__dict__[k] is None}
        assert ppn_none_dict == {}  # 无None的ppn被load出来，不可以有None
        ppn_dict = ppn._to_dict()
        assert ppn_dict.get("pipenode_id") == self.pipenode_id
        del ppn_dict["pipenode_id"]
        assert ppn_dict == self.answer_ppn_dict_partial
        flag, query_data = PipeNodeInfo().query_by_id(self.pipenode_id)
        assert flag is True
        assert isinstance(query_data.get("create_time"), str)
        assert isinstance(query_data.get("last_write_time"), str)
        del query_data["create_time"]
        del query_data["last_write_time"]
        assert query_data.get("pipenode_id") == self.pipenode_id
        del query_data["pipenode_id"]
        assert query_data == self.answer_ppn_dict_partial

    def test_pipenode_create_query_and_delete(self):
        # create
        ppn = PipeNode(conf=self.workflow_conf_dict)
        ppn_id = ppn.ppn_id
        self.appeared_id_list.append(ppn_id)
        file_path = os.path.join(self.pipenode_db_folder, ppn_id + ".pkl")
        assert os.path.exists(file_path)
        # query
        flag, query_data = PipeNodeInfo().query_by_id(ppn_id)
        assert flag is True
        assert isinstance(query_data.get("pipenode_id"), str)
        assert query_data.get("pipenode_id") == ppn_id
        del query_data["pipenode_id"]
        assert isinstance(query_data.get("create_time"), str)
        assert isinstance(query_data.get("last_write_time"), str)
        del query_data["create_time"]
        del query_data["last_write_time"]
        assert query_data == self.answer_ppn_dict_partial
        # delete
        flag, msg = PipeNodeInfo().delete_by_id(ppn_id)
        assert flag is True
        assert msg is True
        assert not os.path.exists(file_path)

    def test_pipenode_outputs_r(self):
        # create
        ppn = PipeNode(conf=self.workflow_conf_dict)
        ppn_id = ppn.ppn_id
        self.appeared_id_list.append(ppn_id)
        file_path = os.path.join(self.pipenode_db_folder, ppn_id + ".pkl")
        assert os.path.exists(file_path)
        # outputs_r赋值 update
        outputs_r = {"results": 1}
        update_dict = {"outputs_r": outputs_r}
        flag, msg = PipeNodeInfo().update_by_id(ppn_id, update_dict)
        assert flag is True
        assert msg is True
        # outputs_r取值
        flag, query_data = PipeNodeInfo().query_by_id(ppn_id)
        assert flag is True
        assert isinstance(query_data.get("pipenode_id"), str)
        assert query_data.get("pipenode_id") == ppn_id
        del query_data["pipenode_id"]
        assert isinstance(query_data.get("create_time"), str)
        assert isinstance(query_data.get("last_write_time"), str)
        del query_data["create_time"]
        del query_data["last_write_time"]
        assert query_data.get("outputs_r") == outputs_r
        # outputs_r读档
        load_ppn = PipeNode().load_by_id(ppn_id)
        assert load_ppn.ppn_id == ppn_id
        assert load_ppn.outputs_r == outputs_r
        # delete
        flag, msg = PipeNodeInfo().delete_by_id(ppn_id)
        assert flag is True
        assert msg is True
        assert not os.path.exists(file_path)
