# -*- coding:utf8 -*-
"""
测试src/utils/utils.py下所有功能是否正确执行
这些都是零散的小函数，如果有需要，随时可以独立放到外面
"""


import pytest
import time
import uuid
from utils.log import get_logger
from utils.utils import PipeTaskStatus, new_id, time_2_asctime, asctime_2_time


logger = get_logger(__name__)


class TestPipeTaskStatus(object):

    def test_status_names(self):
        assert PipeTaskStatus.PREPARATION.name == "PREPARATION"
        assert PipeTaskStatus.DOING.name == "DOING"
        assert PipeTaskStatus.SUCCESS.name == "SUCCESS"
        assert PipeTaskStatus.FAIL.name == "FAIL"
        assert PipeTaskStatus.RESTARTING.name == "RESTARTING"

    def test_status_next(self):
        assert PipeTaskStatus.PREPARATION.value == ["DOING", "FAIL"]
        assert PipeTaskStatus.DOING.value == ["SUCCESS", "FAIL", "RESTARTING"]
        assert PipeTaskStatus.SUCCESS.value == []
        assert PipeTaskStatus.FAIL.value == ["RESTARTING"]
        assert PipeTaskStatus.RESTARTING.value == ["DOING", "FAIL", None]

    def test_status_pre(self):
        eval_name = "PipeTaskStatus.{}.name"
        eval_value = "PipeTaskStatus.{}.value"
        pre_status_dict = {"PREPARATION": [],
                           "DOING": ["PREPARATION", "RESTARTING"],
                           "SUCCESS": ["DOING"],
                           "FAIL": ["PREPARATION", "DOING", "RESTARTING"],
                           "RESTARTING": ["FAIL", "DOING"]}
        for now_status in pre_status_dict:
            for pre_status in pre_status_dict.get(now_status):
                if now_status is not None and pre_status is not None:
                    assert eval(eval_name.format(now_status)) in eval(eval_value.format(pre_status))


class TestId(object):

    def test_new_id(self):
        my_id = new_id(is_log=False)
        uu_id = str(uuid.uuid4())
        assert isinstance(my_id, str)
        assert len(my_id) == len(uu_id)


class TestTime(object):

    def setup_class(self):
        # asctime就是转化不到time.time的精度
        self.test_time = 1645695438.6341782
        self.test_asctime = 'Thu Feb 24 17:37:18 2022'
        self.test_str_time_1 = time.localtime(self.test_time)
        self.test_str_time_2 = time.strptime(self.test_asctime, "%a %b %d %H:%M:%S %Y")
        self.mk_time_1 = 1645695438.0
        self.mk_time_2 = 1645695438.0

    def teardown_class(self):
        pass

    def test_time_2_asctime(self):
        asc_time = time_2_asctime(self.test_time)
        assert asc_time == self.test_asctime

    def test_asctime_2_time(self):
        my_time = asctime_2_time(self.test_asctime)
        assert my_time == self.mk_time_1

    def test_time_more(self):
        assert self.test_str_time_1 != self.test_str_time_2

        my_time = time.time()
        my_time_back = asctime_2_time(time_2_asctime(my_time))
        assert 0 < my_time - my_time_back < 1

        my_asctime = time.asctime()
        my_asctime_back = time_2_asctime(asctime_2_time(my_asctime))
        assert my_asctime == my_asctime_back
