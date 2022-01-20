# -*- coding:utf8 -*-
"""
测试pipeline能否被成功实例化，检查，调用
"""


import os
import pytest
import shutil
from db.typing import PipeLineInfo
from db.typing import PipeNodeInfo
from pipeline.pipeline import PipeLine
from utils.log import get_logger
from utils.io import Path
from utils.json_op import Json


logger = get_logger(__name__)


class TestPipeLine(object):

    def setup_class(self,):
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_dag_workflow_file_path = os.path.join(self.fake_path, "fake_dag_workflow.json")
        self.fake_workflow_conf = Json.file_to_json_without_comments(self.fake_dag_workflow_file_path)
        self.answer_calc_dag = {"f1": {"next_nodes": ["f2"], "prep_nodes": []},
                                "f2": {"next_nodes": ["f3", "f5"], "prep_nodes": ["f1"]},
                                "f3": {"next_nodes": ["f4"], "prep_nodes": ["f2"]},
                                "f4": {"next_nodes": ["f6"], "prep_nodes": ["f3"]},
                                "f5": {"next_nodes": ["f6"], "prep_nodes": ["f2"]},
                                "f6": {"next_nodes": ["f7"], "prep_nodes": ["f4", "f5"]},
                                "f7": {"next_nodes": [], "prep_nodes": ["f6"]}}
        self.answer_topo_order_1 = ["f1", "f2", "f3", "f5", "f4", "f6", "f7"]
        self.answer_topo_order_2 = ["f1", "f2", "f5", "f3", "f4", "f6", "f7"]
        self.appeared_ppl_id_list = []
        self.appeared_ppn_id_list = []
        self.top_path = Path._get_full_path()
        # 先删
        self.pipeline_db_folder = os.path.join(self.top_path, "results", "pipeline_info")
        self.pipenode_db_folder = os.path.join(self.top_path, "results", "pipenode_info")

        if os.path.exists(self.pipenode_db_folder):
            all_file_name_list = os.listdir(self.pipenode_db_folder)
            all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
            for test_file in all_test_ahead_list:
                os.remove(os.path.join(self.pipenode_db_folder, test_file))
        else:
            os.makedirs(self.pipenode_db_folder)

        if os.path.exists(self.pipeline_db_folder):
            all_file_name_list = os.listdir(self.pipeline_db_folder)
            all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
            for test_file in all_test_ahead_list:
                os.remove(os.path.join(self.pipeline_db_folder, test_file))
        else:
            os.makedirs(self.pipeline_db_folder)
        # 再粘
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_path, "test_results", "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_path, "test_results", "fake_pipeline_info")

        all_file_name_list = os.listdir(self.fake_pipenode_pkl_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            fake_pipenode_pkl_path = os.path.join(self.fake_pipenode_pkl_folder, test_file)
            shutil.copy(fake_pipenode_pkl_path, self.pipenode_db_folder)

        all_file_name_list = os.listdir(self.fake_pipeline_pkl_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            fake_pipeline_pkl_path = os.path.join(self.fake_pipeline_pkl_folder, test_file)
            shutil.copy(fake_pipeline_pkl_path, self.pipeline_db_folder)
        self.pipeline_id = "test_load_pipeline_id"

    def teardown_class(self):
        # 删除粘贴过去的文件
        all_file_name_list = os.listdir(self.pipenode_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipenode_db_folder, test_file))
        all_file_name_list = os.listdir(self.pipeline_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipeline_db_folder, test_file))
        # 防止delete失败，尝试删除所有测试中产生的id
        for appeared_id in self.appeared_ppn_id_list:
            file_path = os.path.join(self.pipenode_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                logger.warning("find ppn_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)
        for appeared_id in self.appeared_ppl_id_list:
            file_path = os.path.join(self.pipeline_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                logger.warning("find ppl_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)

    def test_create_blank_and_load_ppl(self):
        self.appeared_ppl_id_list.append(self.pipeline_id)
        ppl = PipeLine().load_by_id(self.pipeline_id)
        ppl_none_dict = {k: ppl.__dict__[k] for k in ppl.__dict__ if ppl.__dict__[k] is None and k != "flags"}
        assert ppl_none_dict == {}  # 无None的ppl被load出来，不可以有None
        assert ppl.ppl_id == self.pipeline_id
        assert ppl.topo_order_list == self.answer_topo_order_1 or ppl.topo_order_list == self.answer_topo_order_2
        assert ppl.dag_dict == self.answer_calc_dag
        flag, query_data = PipeLineInfo().query_by_id(self.pipeline_id)
        assert flag is True
        assert query_data.get("pipeline_id") == self.pipeline_id
        assert isinstance(query_data.get("pipeline_name"), str)
        assert isinstance(query_data.get("create_time"), str)
        assert isinstance(query_data.get("last_write_time"), str)
        assert query_data.get("topo_order_list") == self.answer_topo_order_1 \
            or query_data.get("topo_order_list") == self.answer_topo_order_2
        assert query_data.get("dag_dict") == self.answer_calc_dag

    def test_pipeline_create_query_update_and_delete(self):
        # create
        ppl = PipeLine(workflow_conf=self.fake_workflow_conf)
        ppl_id = ppl.ppl_id
        node_id_dict = ppl.node_id_dict
        self.appeared_ppl_id_list.append(ppl_id)
        self.appeared_ppn_id_list += list(node_id_dict.values())
        assert ppl.topo_order_list == self.answer_topo_order_1 or ppl.topo_order_list == self.answer_topo_order_2
        assert ppl.dag_dict == self.answer_calc_dag
        file_path = os.path.join(self.pipeline_db_folder, ppl_id + ".pkl")
        assert os.path.exists(file_path)
        for ppn_id in node_id_dict.values():
            ppn_file_path = os.path.join(self.pipenode_db_folder, ppn_id + ".pkl")
            os.path.exists(ppn_file_path)
        # query
        flag, query_data = PipeLineInfo().query_by_id(ppl_id)
        assert flag is True
        assert query_data.get("pipeline_id") == ppl_id
        assert isinstance(query_data.get("pipeline_name"), str)
        assert isinstance(query_data.get("create_time"), str)
        assert isinstance(query_data.get("last_write_time"), str)
        assert query_data.get("topo_order_list") == self.answer_topo_order_1 \
            or query_data.get("topo_order_list") == self.answer_topo_order_2
        assert query_data.get("dag_dict") == self.answer_calc_dag
        # delete
        # 注意 pipeline的删除，需要首先删除关联pipenodes！
        for ppn_id in node_id_dict.values():
            flag, msg = PipeNodeInfo().delete_by_id(ppn_id)
            assert flag is True
            assert msg is True
            ppn_file_path = os.path.join(self.pipenode_db_folder, ppn_id + ".pkl")
            assert not os.path.exists(ppn_file_path)
        flag, msg = PipeLineInfo().delete_by_id(ppl_id)
        assert flag is True
        assert msg is True
        assert not os.path.exists(file_path)

    @pytest.mark.skip("pass")
    def test_config_in_pipeline(self):
        pass
