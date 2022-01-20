# -*- coding:utf8 -*-
"""
测试src/db下所有功能是否正确执行
"""


import copy
import numpy as np
import os
import shutil
import time
import pytest
import uuid
from db.local_db import LocalDb
from db.typing import PipeTaskInfo, PipeNodeInfo, PipeLineInfo
from utils.io import Path
from utils.log import get_logger
from utils.utils import PipeTaskStatus


logger = get_logger(__name__)


class TmpTyping(LocalDb):

    def __init__(self):
        super(TmpTyping, self).__init__()
        self.table_type = "tmp_info"
        self.design_table_type.update({  # db会自动添加create_time:str和last_write_time:str两项
            "id": str,

            "fake_str": str,
            "fake_int": int,
            "fake_float": float,
            "fake_bool": bool,
            "fake_list": list,
            "fake_dict": dict,
            "fake_tuple": tuple,
            "fake_object": object,
            "fake_np_int": np.int,
            "fake_np_float": np.float,
            "fake_none": None
        })
        self.map_id = "id"
        self._init_db_folder()


class TestLocalDb(object):
    """
    test class LocalDb functions
    """

    def setup_class(self):
        self.top_path = Path._get_full_path(relative_path="", base_path_type="top")  # 此处没使用config避免循环引用
        self.db_folder = os.path.join(self.top_path, "results", "tmp_info")
        if os.path.exists(self.db_folder):
            shutil.rmtree(self.db_folder)
        self.test_data = {
            "id": str(uuid.uuid4()),

            "fake_str": "str",
            "fake_int": 1,
            "fake_float": 0.01,
            "fake_bool": True,
            "fake_list": [1, 2],
            "fake_dict": {1: 1},
            "fake_tuple": (1,),
            "fake_np_int": np.int(1),
            "fake_np_float": np.float(0.01),
            "fake_none": {"err_msg": "fake"}
        }
        self.test_data_same = copy.deepcopy(self.test_data)
        self.update_data = {"fake_bool": False}
        self.update_data_same = copy.deepcopy(self.update_data)
        self.update_data_wrong = {"fake_bool": False, "id": str(uuid.uuid4())}
        self.test_condition_true = {"id": self.test_data.get("id"), "fake_bool": True}
        self.test_condition_true_same = copy.deepcopy(self.test_condition_true)
        self.test_condition_false = {"id": self.test_data.get("id"), "fake_bool": False}
        self.test_condition_false_same = copy.deepcopy(self.test_condition_false)
        self.test_wrong_data = {
            "id": str(uuid.uuid4()),

            "fake_str": 1
        }
        self.test_wrong_data_same = copy.deepcopy(self.test_wrong_data)
        self.test_wrong_condition = {}

    def teardown_class(self):
        if os.path.exists(self.db_folder):
            shutil.rmtree(self.db_folder)

    def test_001_db_folder_init(self):  # 001使得执行顺序靠前
        assert not os.path.exists(self.db_folder)
        db_info = TmpTyping()
        assert os.path.exists(self.db_folder)

    def test_insert_update_query_and_delete(self):
        db_info = TmpTyping()
        # test insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        # test query
        flag, data = db_info.query(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert data is not False
        for key in self.test_data:
            assert key in data
            assert self.test_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        flag, data_same = db_info.query(self.test_condition_true)
        assert flag
        assert data_same is not False
        assert data == data_same
        assert (data is not data_same)  # 保证它们俩的地址不是一个
        # test update
        time.sleep(2)  # 睡两秒，防止秒级时间update和insert一致
        flag, msg = db_info.update(self.test_condition_false, self.update_data)  # 应该找不到
        logger.info("Warning is supposed here!")
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is False
        flag, msg = db_info.update(self.test_condition_true, self.update_data_wrong)  # 应该找得到但禁止更新
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(msg, str)
        flag, msg = db_info.update(self.test_condition_true, self.update_data)  # 应该得到
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is True
        flag, data_1 = db_info.query(self.test_condition_true)  # 这个就应该是找不到的
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert data_1 is False
        flag, data_1 = db_info.query(self.test_condition_false)
        assert self.test_condition_false == self.test_condition_false_same
        assert flag
        assert data_1 is not False
        assert isinstance(data_1.get("create_time"), str)
        assert isinstance(data_1.get("last_write_time"), str)
        assert data.get("create_time") == data_1.get("create_time")
        assert data.get("last_write_time") != data_1.get("last_write_time")
        assert data.get("fake_bool")
        assert data_1.get("fake_bool") is False
        del data_1["fake_bool"]
        test_data_copy = copy.deepcopy(self.test_data)
        del test_data_copy["fake_bool"]
        for key in test_data_copy:
            assert key in data_1
            assert test_data_copy[key] == data_1.get(key)
        # test delete
        flag, msg = db_info.delete(self.test_condition_true)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is False
        flag, msg = db_info.delete(self.test_condition_false)  # 可以找到
        assert self.test_condition_false == self.test_condition_false_same
        assert flag
        assert msg is True

    def test_insert_wrong(self):
        db_info = TmpTyping()
        test_id = self.test_wrong_data.get("id")
        pkl_path = os.path.join(self.db_folder, test_id + ".pkl")
        flag, msg = db_info.insert(self.test_wrong_data)
        logger.info("Error is supposed here!")
        assert self.test_wrong_data == self.test_wrong_data_same  # 要保证数据不会被变动
        assert not flag
        assert isinstance(msg, str)
        assert not os.path.exists(pkl_path)

    def test_update_wrong(self):
        db_info = TmpTyping()
        # insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        # test update wrong
        wrong_update_data = {
            "id": self.test_data.get("id"),

            "fake_int": "str"
        }
        time.sleep(2)  # 睡两秒，防止秒级时间update和insert一致
        flag, msg = db_info.update(self.test_condition_true, wrong_update_data)  # 应该找不到
        logger.info("Error is supposed here!")
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert not flag
        assert isinstance(msg, str)
        # 再次query，结果应该与初始insert一致，即保证错误的update不会更改原有db
        flag, data = db_info.query(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert data is not False
        for key in self.test_data:
            assert key in data
            assert self.test_data[key] == data.get(key)
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        # delete
        flag, msg = db_info.delete(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True

    def test_query_wrong(self):
        db_info = TmpTyping()
        # insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        # wrong query
        flag, data = db_info.query(self.test_wrong_condition)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(data, str)
        # delete
        flag, msg = db_info.delete(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True

    def test_delete_wrong(self):
        db_info = TmpTyping()
        # insert
        test_id = self.test_data.get("id")
        pkl_path = os.path.join(self.db_folder, test_id + ".pkl")
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        # wrong delete
        flag, msg = db_info.delete(self.test_wrong_condition)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(msg, str)
        assert os.path.exists(pkl_path)
        # delete
        flag, msg = db_info.delete(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True
        assert not os.path.exists(pkl_path)


class TestTyping(object):
    """
    test python file typing.py functions
    all test ids in this case must be test_<id>
    """

    def setup_class(self):
        self.fake_path = Path._get_full_path(relative_path="", base_path_type="top")  # 此处没使用config避免循环引用
        self.db_folder = os.path.join(self.fake_path, "results", "{}")
        self.db_folder_list = [self.db_folder.format(i) for i in ["pipeline_info", "pipenode_info", "pipetask_info"]]
        for db_folder in self.db_folder_list:
            if os.path.exists(db_folder):
                all_file_name_list = os.listdir(db_folder)
                all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
                for test_file in all_test_ahead_list:
                    os.remove(os.path.join(db_folder, test_file))

        self.test_ppt_data = {
            "pipetask_id": "test_" + str(uuid.uuid4()),

            "pipeline_id": str(uuid.uuid4()),
            "finish_node_list": [],
            "pipetask_status": PipeTaskStatus.PREPARATION.name,
            "flags": None
        }
        self.test_ppt_data_same = copy.deepcopy(self.test_ppt_data)
        self.test_ppt_condition = {"pipetask_id": self.test_ppt_data.get("pipetask_id")}
        self.test_ppt_update_data = {"finish_node_list": ["1"]}

        self.test_ppl_data = {
            "pipeline_id": "test_" + str(uuid.uuid4()),
            "pipeline_name": "test_pipeline",

            "dag_dict": {},
            "topo_order_list": [],
            "config": None,
            "node_id_dict": {},
            "flags": None
        }
        self.test_ppl_data_same = copy.deepcopy(self.test_ppl_data)
        self.test_ppl_condition = {"pipeline_id": self.test_ppl_data.get("pipeline_id")}
        self.test_ppl_update_data = {"topo_order_list": ["1"]}

        self.test_ppn_data = {
            "pipenode_id": "test_" + str(uuid.uuid4()),
            "pipenode_name": "test_pipenode",

            "func_des": ["None.fake.fake_core", "fake_function", ""],
            "func_str": "fake",
            "type": "cold",
            "inputs": [],
            "outputs": [],
            "next_nodes": [],
            "prep_nodes": [],
            "outputs_r": {},
            "flags": ""
        }
        self.test_ppn_data_same = copy.deepcopy(self.test_ppn_data)
        self.test_ppn_condition = {"pipenode_id": self.test_ppn_data.get("pipenode_id")}
        self.test_ppn_update_data = {"inputs": ["1"]}

    def teardown_class(self):
        for db_folder in self.db_folder_list:
            if os.path.exists(db_folder):
                all_file_name_list = os.listdir(db_folder)
                all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
                for test_file in all_test_ahead_list:
                    os.remove(os.path.join(db_folder, test_file))

    def test_pipetaskinfo(self):
        db_info = PipeTaskInfo()
        # test insert
        flag, msg = db_info.insert(self.test_ppt_data)
        assert flag
        assert msg is True
        file_path = os.path.join(self.db_folder.format("pipetask_info"), self.test_ppt_data.get("pipetask_id") + ".pkl")
        assert os.path.exists(file_path)
        # test query
        flag, data = db_info.query(self.test_ppt_condition)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.test_ppt_data_same
        assert os.path.exists(file_path)
        # test update
        time.sleep(2)
        flag, msg = db_info.update(self.test_ppt_condition, self.test_ppt_update_data)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query(self.test_ppt_condition)
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("finish_node_list") == self.test_ppt_update_data.get("finish_node_list")
        # test delete
        flag, msg = db_info.delete(self.test_ppt_condition)
        assert flag
        assert msg is True
        assert not os.path.exists(file_path)

    def test_pipenodeinfo(self):
        db_info = PipeNodeInfo()
        # test insert
        flag, msg = db_info.insert(self.test_ppn_data)
        assert flag
        assert msg is True
        file_path = os.path.join(self.db_folder.format("pipenode_info"), self.test_ppn_data.get("pipenode_id") + ".pkl")
        assert os.path.exists(file_path)
        # test query
        flag, data = db_info.query(self.test_ppn_condition)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.test_ppn_data_same
        assert os.path.exists(file_path)
        # test update
        time.sleep(2)
        flag, msg = db_info.update(self.test_ppn_condition, self.test_ppn_update_data)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query(self.test_ppn_condition)
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("inputs") == self.test_ppn_update_data.get("inputs")
        # test delete
        flag, msg = db_info.delete(self.test_ppn_condition)
        assert flag
        assert msg is True
        assert not os.path.exists(file_path)

    def test_pipelineinfo(self):
        db_info = PipeLineInfo()
        # test insert
        flag, msg = db_info.insert(self.test_ppl_data)
        assert flag
        assert msg is True
        file_path = os.path.join(self.db_folder.format("pipeline_info"), self.test_ppl_data.get("pipeline_id") + ".pkl")
        assert os.path.exists(file_path)
        # test query
        flag, data = db_info.query(self.test_ppl_condition)
        assert flag
        assert isinstance(data.get("create_time"), str)
        assert isinstance(data.get("last_write_time"), str)
        create_time_1 = data.get("create_time")
        last_write_time_1 = data.get("last_write_time")
        del data["create_time"]
        del data["last_write_time"]
        assert data == self.test_ppl_data_same
        assert os.path.exists(file_path)
        # test update
        time.sleep(2)
        flag, msg = db_info.update(self.test_ppl_condition, self.test_ppl_update_data)
        assert flag
        assert msg is True
        flag, data_2 = db_info.query(self.test_ppl_condition)
        assert data_2.get("create_time") == create_time_1
        assert data_2.get("last_write_time") != last_write_time_1
        assert data_2.get("topo_order_list") == self.test_ppl_update_data.get("topo_order_list")
        # test delete
        flag, msg = db_info.delete(self.test_ppl_condition)
        assert flag
        assert msg is True
        assert not os.path.exists(file_path)
