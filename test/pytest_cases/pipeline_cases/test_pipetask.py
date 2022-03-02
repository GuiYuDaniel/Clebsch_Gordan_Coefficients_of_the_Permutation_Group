# -*- coding:utf8 -*-
"""
测试pipetask能否顺利实例化，执行
"""


import os
import pytest
import shutil
from db.typing import PipeNodeInfo, PipeTaskInfo
from pipeline.pipeline import PipeLine
from pipeline.pipetask import PipeTask
from utils.log import get_logger
from utils.io import Path
from utils.json_op import Json
from utils.utils import PipeTaskStatus


logger = get_logger(__name__)


# @pytest.mark.skip("pass")
class TestPipeTask(object):
    """
    主要是测试pipetask的流程问题，使用一个具有代表性，尽可能全面的workflow来测
    start，restart都要测
    """

    def setup_class(self):
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_dag_workflow_file_path = os.path.join(self.fake_path, "fake_dag_workflow.json")
        self.fake_workflow_conf = Json.file_to_json_without_comments(self.fake_dag_workflow_file_path)
        self.finish_node_list_before_start = []
        self.answer_topo_order_1 = ["f1", "f2", "f3", "f5", "f4", "f6", "f7"]
        self.answer_topo_order_2 = ["f1", "f2", "f5", "f3", "f4", "f6", "f7"]
        self.first_input_args_before_start = None
        self.first_input_args_after_start = tuple([])
        self.first_input_kwargs_before_start = None
        self.first_input_kwargs_after_start = {"flag": True}
        self.ppt_status_before_start = PipeTaskStatus.PREPARATION.name
        self.ppt_status_after_start = PipeTaskStatus.SUCCESS.name
        self.appeared_ppl_id_list = []
        self.appeared_ppn_id_list = []
        self.appeared_ppt_id_list = []
        self.top_path = Path._get_full_path()
        # 先删
        self.pipeline_db_folder = os.path.join(self.top_path, "results", "pipeline_info")
        self.pipenode_db_folder = os.path.join(self.top_path, "results", "pipenode_info")
        self.pipetask_db_folder = os.path.join(self.top_path, "results", "pipetask_info")

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

        if os.path.exists(self.pipetask_db_folder):
            all_file_name_list = os.listdir(self.pipetask_db_folder)
            all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
            for test_file in all_test_ahead_list:
                os.remove(os.path.join(self.pipetask_db_folder, test_file))
        else:
            os.makedirs(self.pipetask_db_folder)
        # 再粘
        self.fake_pipenode_pkl_folder = os.path.join(self.fake_path, "test_results", "fake_pipenode_info")
        self.fake_pipeline_pkl_folder = os.path.join(self.fake_path, "test_results", "fake_pipeline_info")
        self.fake_pipetask_pkl_folder = os.path.join(self.fake_path, "test_results", "fake_pipetask_info")
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
        all_file_name_list = os.listdir(self.fake_pipetask_pkl_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            fake_pipetask_pkl_path = os.path.join(self.fake_pipetask_pkl_folder, test_file)
            shutil.copy(fake_pipetask_pkl_path, self.pipetask_db_folder)
        self.test_ppt_restart_ppn_id_list = [
        ]
        self.test_ppt_restart_id_1 = "test_pipetask_restart_1_id"
        self.test_ppt_restart_id_2 = "test_pipetask_restart_2_id"
        self.test_ppt_restart_id_3 = "test_pipetask_restart_3_id"

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
        all_file_name_list = os.listdir(self.pipetask_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipetask_db_folder, test_file))
        # 防止delete失败，尝试删除所有测试中产生的id
        for appeared_id in self.appeared_ppn_id_list:
            file_path = os.path.join(self.pipenode_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppn_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)
        for appeared_id in self.appeared_ppl_id_list:
            file_path = os.path.join(self.pipeline_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppl_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)
        for appeared_id in self.appeared_ppt_id_list:
            file_path = os.path.join(self.pipetask_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppt_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)

    def test_pipetask_create(self):
        ppl = PipeLine(workflow_conf=self.fake_workflow_conf,
                       ppl_name="test_pipeline_all",
                       ppl_id="test_pipeline_all_id")
        ppl_id = ppl.ppl_id
        node_id_dict = ppl.node_id_dict
        self.appeared_ppl_id_list.append(ppl_id)
        self.appeared_ppn_id_list += list(node_id_dict.values())
        ppt = PipeTask(ppl=ppl, use_name_replace_id="test_pipeline_all_id")
        ppt_id = ppt.ppt_id
        self.appeared_ppt_id_list.append(ppt_id)
        assert ppl_id == ppt.ppl_id
        assert ppt.finish_node_list == self.finish_node_list_before_start
        assert ppt.first_input_args == self.first_input_args_before_start
        assert ppt.first_input_kwargs == self.first_input_kwargs_before_start
        assert ppt.ppt_status == self.ppt_status_before_start
        ppt.start(flag=True)
        assert ppt.finish_node_list == self.answer_topo_order_1 \
            or ppt.finish_node_list == self.answer_topo_order_2
        assert ppt.first_input_args == self.first_input_args_after_start
        assert ppt.first_input_kwargs == self.first_input_kwargs_after_start
        assert ppt.ppt_status == self.ppt_status_after_start
        file_path = os.path.join(self.pipetask_db_folder, ppt_id + ".pkl")
        assert os.path.exists(file_path)
        flag, query_data = PipeTaskInfo().query_by_id(ppt_id)
        assert flag is True
        assert query_data.get("pipetask_id") == ppt_id
        assert query_data.get("pipeline_id") == ppl_id
        assert query_data.get("finish_node_list") == ppt.finish_node_list

    def test_pipetask_restart_1(self):
        """stop at first node(1)
        [1][2][354][6][7]"""
        ppt = PipeTask().load_by_id(self.test_ppt_restart_id_1)
        ppt_none_dict = {k: ppt.__dict__[k] for k in ppt.__dict__ if ppt.__dict__[k] is None and k != "flags"}
        assert ppt_none_dict == {}  # 无None的ppt被load出来，不可以有None
        ppt_id = ppt.ppt_id
        ppl_id = ppt.ppl_id
        self.appeared_ppt_id_list.append(ppt_id)
        self.appeared_ppl_id_list.append(ppl_id)
        ppt.restart()
        assert ppt.finish_node_list == self.answer_topo_order_1 \
               or ppt.finish_node_list == self.answer_topo_order_2
        assert ppt.first_input_args == self.first_input_args_after_start
        assert ppt.first_input_kwargs == self.first_input_kwargs_after_start
        assert ppt.ppt_status == self.ppt_status_after_start
        file_path = os.path.join(self.pipetask_db_folder, ppt_id + ".pkl")
        assert os.path.exists(file_path)
        flag, query_data = PipeTaskInfo().query_by_id(ppt_id)
        assert flag is True
        assert query_data.get("pipetask_id") == ppt_id
        assert query_data.get("pipeline_id") == ppl_id
        assert query_data.get("finish_node_list") == ppt.finish_node_list

    def test_pipetask_restart_2(self):
        """stop at middle node(5)
        [1][2][354][6][7]"""
        ppt = PipeTask().load_by_id(self.test_ppt_restart_id_2)
        ppt_none_dict = {k: ppt.__dict__[k] for k in ppt.__dict__ if ppt.__dict__[k] is None and k != "flags"}
        assert ppt_none_dict == {}  # 无None的ppt被load出来，不可以有None
        ppt_id = ppt.ppt_id
        ppl_id = ppt.ppl_id
        self.appeared_ppt_id_list.append(ppt_id)
        self.appeared_ppl_id_list.append(ppl_id)
        ppt.restart()
        assert ppt.finish_node_list == self.answer_topo_order_1 \
               or ppt.finish_node_list == self.answer_topo_order_2
        assert ppt.first_input_args == self.first_input_args_after_start
        assert ppt.first_input_kwargs == self.first_input_kwargs_after_start
        assert ppt.ppt_status == self.ppt_status_after_start
        file_path = os.path.join(self.pipetask_db_folder, ppt_id + ".pkl")
        assert os.path.exists(file_path)
        flag, query_data = PipeTaskInfo().query_by_id(ppt_id)
        assert flag is True
        assert query_data.get("pipetask_id") == ppt_id
        assert query_data.get("pipeline_id") == ppl_id
        assert query_data.get("finish_node_list") == ppt.finish_node_list

    def test_pipetask_restart_3(self):
        """stop at last node(7)
        [1][2][354][6][7]"""
        ppt = PipeTask().load_by_id(self.test_ppt_restart_id_3)
        ppt_none_dict = {k: ppt.__dict__[k] for k in ppt.__dict__ if ppt.__dict__[k] is None and k != "flags"}
        assert ppt_none_dict == {}  # 无None的ppt被load出来，不可以有None
        ppt_id = ppt.ppt_id
        ppl_id = ppt.ppl_id
        self.appeared_ppt_id_list.append(ppt_id)
        self.appeared_ppl_id_list.append(ppl_id)
        ppt.restart()
        assert ppt.finish_node_list == self.answer_topo_order_1 \
               or ppt.finish_node_list == self.answer_topo_order_2
        assert ppt.first_input_args == self.first_input_args_after_start
        assert ppt.first_input_kwargs == self.first_input_kwargs_after_start
        assert ppt.ppt_status == self.ppt_status_after_start
        file_path = os.path.join(self.pipetask_db_folder, ppt_id + ".pkl")
        assert os.path.exists(file_path)
        flag, query_data = PipeTaskInfo().query_by_id(ppt_id)
        assert flag is True
        assert query_data.get("pipetask_id") == ppt_id
        assert query_data.get("pipeline_id") == ppl_id
        assert query_data.get("finish_node_list") == ppt.finish_node_list


# @pytest.mark.skip("pass")
class TestSpecialWorkflow(object):
    """
    一些特别的function，写成workflow，在这里测试
    """

    def setup_class(self):
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.fake_wfl_output_str = os.path.join(self.fake_path, "fake_workflow_output_str.json")
        self.fake_wfl_output_str_conf = Json.file_to_json_without_comments(self.fake_wfl_output_str)
        # self.fake_wfl_output_str = os.path.join(self.fake_path, "fake_workflow_output_str.json")
        # self.fake_wfl_output_str_conf = Json.file_to_json_without_comments(self.fake_wfl_output_str)
        self.appeared_ppl_id_list = []
        self.appeared_ppn_id_list = []
        self.appeared_ppt_id_list = []
        self.top_path = Path._get_full_path()
        self.pipeline_db_folder = os.path.join(self.top_path, "results", "pipeline_info")
        self.pipenode_db_folder = os.path.join(self.top_path, "results", "pipenode_info")
        self.pipetask_db_folder = os.path.join(self.top_path, "results", "pipetask_info")
        all_file_name_list = os.listdir(self.pipenode_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipenode_db_folder, test_file))
        all_file_name_list = os.listdir(self.pipeline_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipeline_db_folder, test_file))
        all_file_name_list = os.listdir(self.pipetask_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipetask_db_folder, test_file))

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
        all_file_name_list = os.listdir(self.pipetask_db_folder)
        all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        for test_file in all_test_ahead_list:
            os.remove(os.path.join(self.pipetask_db_folder, test_file))
        # 防止delete失败，尝试删除所有测试中产生的id
        for appeared_id in self.appeared_ppn_id_list:
            file_path = os.path.join(self.pipenode_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppn_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)
        for appeared_id in self.appeared_ppl_id_list:
            file_path = os.path.join(self.pipeline_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppl_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)
        for appeared_id in self.appeared_ppt_id_list:
            file_path = os.path.join(self.pipetask_db_folder, appeared_id + ".pkl")
            if os.path.exists(file_path):
                # logger.warning("find ppt_id={} still exists, remove it".format(appeared_id))
                os.remove(file_path)

    def test_output_str_func(self):
        ppl = PipeLine(workflow_conf=self.fake_wfl_output_str_conf, ppl_id="test_output_str_id")
        ppl_id = ppl.ppl_id
        node_id_dict = ppl.node_id_dict
        topo_order_list = ppl.topo_order_list
        self.appeared_ppl_id_list.append(ppl_id)
        self.appeared_ppn_id_list += list(node_id_dict.values())
        ppt = PipeTask(ppl=ppl, use_name_replace_id="test_output_str_id")
        ppt_id = ppt.ppt_id
        self.appeared_ppt_id_list.append(ppt_id)
        ppt.start()  # 对于首节点无参数输入的情况，可以不写参数，首节点conf里的inputs，可以写[]，或[None]，反正没人看
        ppn_1_id = node_id_dict.get(topo_order_list[0])
        ppn_2_id = node_id_dict.get(topo_order_list[1])
        flag, query_data = PipeNodeInfo().query_by_id(ppn_1_id)
        assert flag is True
        assert query_data.get("outputs_r").get("string") == "this is str"
        flag, query_data = PipeNodeInfo().query_by_id(ppn_2_id)
        assert flag is True
        assert query_data.get("outputs_r").get("param1:::flag") is True


# @pytest.mark.skip("pass")
class TestFunc(object):
    """
    对于pipetask可能执行的function，做一下test，保证行为要与预期一致，
    这里仅抽出核心机制测
    """

    def setup_class(self):
        pass
        # self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        # self.fake_wfl_output_str = os.path.join(self.fake_path, "fake_workflow_output_str.json")
        # self.fake_wfl_output_str_conf = Json.file_to_json_without_comments(self.fake_wfl_output_str)
        # self.fake_wfl_output_str = os.path.join(self.fake_path, "fake_workflow_output_str.json")
        # self.fake_wfl_output_str_conf = Json.file_to_json_without_comments(self.fake_wfl_output_str)

    def test_easy_output_str_func_1(self):
        # 这个是抽出来的核心机制
        exec("from test.fake.fake_core import output_str_func")
        exec("from test.fake.fake_core import is_input_str_func")
        func_1 = eval("output_str_func")
        func_2 = eval("is_input_str_func")
        output_1 = func_1()
        output_1_tuple = tuple([output_1])
        output_1_dict = {}
        output_2 = func_2(*output_1_tuple, **output_1_dict)
        assert output_2

    def test_easy_output_str_func_2(self):
        # 这个是抽出来的核心机制
        exec("from test.fake.fake_core import output_str_func")
        exec("from test.fake.fake_core import is_input_str_func")
        func_1 = eval("output_str_func")
        func_2 = eval("is_input_str_func")
        output_1 = func_1()
        output_1_tuple = tuple([])
        output_1_dict = {"string": output_1}
        output_2 = func_2(*output_1_tuple, **output_1_dict)
        assert output_2

    def test_easy_output_all_str_fun(self):
        # 这个是抽出来的核心机制
        exec("from test.fake.fake_core import output_all_str_func")
        exec("from test.fake.fake_core import is_input_all_str_func")
        func_1 = eval("output_all_str_func")
        func_2 = eval("is_input_all_str_func")
        output_1 = func_1(None)
        output_1 = tuple([output_1[0], output_1[1]])
        output_2 = func_2(*output_1, **{})
        assert output_2

    # def test_output_number_func(self):
    #     pass
    #
    # def test_output_list_func(self):
    #     pass
    #
    # def test_output_dict_func(self):
    #     pass
    #
    # def test_output_tuple_func(self):
    #     # python的多变量返回也是tuple，避免和单变量就是tuple情况混淆
    #     pass

    def test_output_multi_str_func_1(self):
        def output_str_func():
            return "this is str 1", "this is str 2"

        def input_str_func(a, b):
            return True if isinstance(a, str) and isinstance(b, str) else False

        func_1 = output_str_func
        func_2 = input_str_func
        output_1 = func_1()
        output_2 = func_2(*output_1)
        assert output_2

    def test_output_multi_str_fun_2(self):
        def output_str_func():
            return "this is str 1", "this is str 2"

        def input_str_func(*args):
            return True if isinstance(args[0], str) and isinstance(args[1], str) else False

        func_1 = output_str_func
        func_2 = input_str_func
        output_1 = func_1()
        output_2 = func_2(*output_1)
        assert output_2

    # def test_output_multi_number_func(self):
    #     pass
    #
    # def test_output_multi_list_func(self):
    #     pass
    #
    # def test_output_multi_dict_func(self):
    #     pass
    #
    # def test_output_multi_tuple_func(self):
    #     pass
    #
    # def test_output_multi_type_func(self):
    #     pass
