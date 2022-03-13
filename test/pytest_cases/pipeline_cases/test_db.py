# -*- coding:utf8 -*-
"""
测试src/db下所有功能是否正确执行
"""


import copy
import numpy as np
import os
import time
import pytest
import uuid
from db.local_db import LocalDb
from db.local_db_protector import DBProtector
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

    def setup(self):
        self.protector = DBProtector("results", extension_name=".test_db_protected")
        self.protector.protector_setup()

        self.top_path = Path._get_full_path(relative_path="", base_path_type="top")  # 此处没使用config避免循环引用
        self.db_folder = os.path.join(self.top_path, "results", "tmp_info")
        self.test_data = {
            "id": str(uuid.uuid4()),

            "fake_str": "str",
            "fake_int": 1,
            "fake_float": 0.01,
            "fake_bool": True,
            "fake_list": [1, 2],
            "fake_dict": {1: 1, "2": "2"},
            "fake_tuple": (1,),
            "fake_np_int": np.int(1),
            "fake_np_float": np.float(0.01),
            "fake_none": {"err_msg": "fake"}
        }
        self.test_data_id = self.test_data.get("id")
        self.deep_partial_dict = {"fake_dict": {1: 1, "3": "3"}}
        self.deep_partial_dict_same = copy.deepcopy(self.deep_partial_dict)
        self.deep_updated_dict = {1: 1, "3": "3"}
        self.deep_not_updated_dict = {1: 1, "2": "2", "3": "3"}
        self.file_name_pkl = os.path.join(self.db_folder, self.test_data.get("id") + ".pkl")
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
        self.test_wrong_condition_1 = []

    def teardown(self):
        self.protector.protector_teardown()

    def test_001_db_folder_init(self):  # 001使得执行顺序靠前
        assert not os.path.exists(self.db_folder)
        db_info = TmpTyping()
        assert os.path.exists(self.db_folder)

    def test_insert_update_query_and_delete(self):
        db_info = TmpTyping()
        assert not os.path.exists(self.file_name_pkl)

        # test insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)

        # insert again
        m_time_before = os.stat(self.file_name_pkl).st_mtime
        flag, msg = db_info.insert(self.test_data)
        logger.info("Error is supposed here!")
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag is False
        assert msg
        assert os.path.exists(self.file_name_pkl)
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before == m_time_after

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

        flag, data = db_info.query_by_id(self.test_data_id)
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
        time.sleep(1.1)  # 睡两秒，防止秒级时间update和insert一致

        flag, msg = db_info.update(self.test_condition_false, self.update_data)  # 应该找不到
        logger.info("Warning is supposed here!")
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert self.test_condition_false == self.test_condition_false_same
        assert flag
        assert msg is False

        flag, msg = db_info.update(self.test_condition_true, self.update_data_wrong)  # 应该找得到但禁止更新
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(msg, str)

        m_time_before = os.stat(self.file_name_pkl).st_mtime
        flag, msg = db_info.update(self.test_condition_true, self.update_data)  # 应该得到
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before != m_time_after

        m_time_before = os.stat(self.file_name_pkl).st_mtime  # update again
        flag, msg = db_info.update_by_id(self.test_data_id, self.update_data)  # 应该得到
        assert self.update_data == self.update_data_same
        assert self.test_data == self.test_data_same
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before == m_time_after

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
        assert os.path.exists(self.file_name_pkl)
        flag, msg = db_info.delete(self.test_condition_true)  # 应该是找不到的
        logger.info("Warning is supposed here!")
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)

        flag, msg = db_info.delete(self.test_condition_false)  # 可以找到
        assert self.test_condition_false == self.test_condition_false_same
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)

        flag, msg = db_info.delete(self.test_condition_false)  # delete again
        logger.info("Warning is supposed here!")
        assert self.test_condition_false == self.test_condition_false_same
        assert flag
        assert msg is False
        assert not os.path.exists(self.file_name_pkl)

    def test_update_deep(self):
        """
        对于python，
        a = {1: 2}
        a.update({3: 4}) == None
        a == {1: 2, 3: 4}

        b = {"data": {1: 2}}
        b.update({"data": {3: 4}}) == None
        b == {"data": {3: 4}}
        不是 不是 不是 {"data": {1: 2, 3: 4}} 哦！！！
        也就是说，当前的localdb，只能update到第一层key-value，不支持深层update
        """
        db_info = TmpTyping()
        assert not os.path.exists(self.file_name_pkl)

        # insert
        flag, msg = db_info.insert(self.test_data)
        assert self.test_data == self.test_data_same  # 要保证数据不会被变动
        assert flag
        assert msg is True
        assert os.path.exists(self.file_name_pkl)

        flag, data = db_info.query(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert isinstance(data, dict)
        assert data.get("fake_dict") == self.test_data_same.get("fake_dict")

        # update deep
        m_time_before = os.stat(self.file_name_pkl).st_mtime
        flag, msg = db_info.update_by_id(self.test_data_id, self.deep_partial_dict)  # 应该得到
        assert self.test_data == self.test_data_same
        assert self.test_condition_true == self.test_condition_true_same
        assert self.deep_partial_dict == self.deep_partial_dict_same
        assert flag
        assert msg is True
        m_time_after = os.stat(self.file_name_pkl).st_mtime
        assert m_time_before != m_time_after
        assert os.path.exists(self.file_name_pkl)

        flag, data_1 = db_info.query(self.test_condition_true)
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert isinstance(data_1, dict)
        assert data_1.get("fake_dict") == self.deep_updated_dict

        # test delete
        assert os.path.exists(self.file_name_pkl)
        flag, msg = db_info.delete(self.test_data)  # 应该是找得到，但不符合条件，因为test_data已经被update了
        logger.info("Warning is supposed here!")
        assert self.test_data == self.test_data_same
        assert flag
        assert msg is False
        assert os.path.exists(self.file_name_pkl)

        flag, msg = db_info.delete_by_id(self.test_data_id)  # 可以找到
        assert self.test_condition_true == self.test_condition_true_same
        assert flag
        assert msg is True
        assert not os.path.exists(self.file_name_pkl)

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
        # time.sleep(1.1)  # 睡1.1秒，防止秒级时间update和insert一致
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

        flag, data = db_info.query(self.test_wrong_condition_1)
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
    TestLocalDB中已经测过基础操作和基础报错了。所以，这里只测正向使用即可
    """

    def setup_class(self):
        self.protector = DBProtector("results", extension_name=".test_db_protected")
        self.protector.protector_setup()

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
        # for db_folder in self.db_folder_list:
        #     if os.path.exists(db_folder):
        #         all_file_name_list = os.listdir(db_folder)
        #         all_test_ahead_list = [i for i in all_file_name_list if i.startswith("test_")]
        #         for test_file in all_test_ahead_list:
        #             os.remove(os.path.join(db_folder, test_file))

        self.protector.protector_teardown()

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
        time.sleep(1.1)
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
        time.sleep(1.1)
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
        time.sleep(1.1)
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
