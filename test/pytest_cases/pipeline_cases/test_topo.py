# -*- coding:utf8 -*-
"""
测试src/utils/topo功能是否正确执行
"""


import os
from utils.json_op import Json
from utils.io import Path
from utils.topo import calc_dag, calc_topo_order


class TestTopo(object):
    """
    test topo functions
    """

    def setup_class(self):
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

    def test_calc_dag(self):
        answer = self.answer_calc_dag
        dag_dict = calc_dag(self.fake_workflow_conf)
        assert dag_dict == answer

    def test_calc_topo_order(self):
        # 1 2 都是可能且合理的结果
        answer_1 = self.answer_topo_order_1
        answer_2 = self.answer_topo_order_2
        topo_order = calc_topo_order(self.answer_calc_dag)
        assert (topo_order == answer_1 or topo_order == answer_2)
