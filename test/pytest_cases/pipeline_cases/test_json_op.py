# -*- coding:utf8 -*-
"""
测试src/utils/json_op功能是否正确执行
"""


import os
import pytest
from utils.io import Path
from utils.json_op import Json


class TestJson(object):
    """
    test class Json functions
    """

    def setup_class(self):
        self.fake_path = Path._get_full_path(relative_path="fake", base_path_type="test")  # 此处没使用config避免循环引用
        self.whereami_file_path = os.path.join(self.fake_path, "fake_workflow_whereami.json")
        self.answer_whereami = [{"name": "first_node",
                                 "func": ["None.test.fake.fake_core", "where_am_i", ""],
                                 "type": "cold",
                                 "inputs": ["point_path"],
                                 "outputs": ["results"],
                                 "extra_args": [],
                                 "extra_kwargs": {},
                                 "next_nodes": [],
                                 "flags": []}]
        self.answer_workflow1_file_path = os.path.join(self.fake_path, "fake_json_with_comments1.json")
        self.answer_workflow1 = \
            [{"name": "first_node",
              "func": ["None.test.fake.fake_core", "where_am_i", ""],
              "type": None,
              "inputs": ["input:cmd_params"],
              "outputs": ["output:results"],
              "extra_args": [],
              "extra_kwargs": {},
              "flags": [],
              "keep1": "### this should keep ###",
              "keep2": ["this should keep"],
              "fix1": {
                       "keep3": "### this should keep ###",
                       "keep4": ["this should be keep",

                                 {}]}},
             ["this should be keep"],
             ]
        self.answer_workflow1_with_comments = \
            ["### a fix dict start 1 ###",
             {"name": "first_node",
              "func": ["None.test.fake.fake_core", "where_am_i", ""],
              "type": None,
              "inputs": ["input:cmd_params"],
              "outputs": ["output:results"],
              "extra_args": [],
              "extra_kwargs": {},
              "flags": [],
              "### delete1 ###": "this should be delete after without comments",
              "keep1": "### this should keep ###",
              "### delete2 ###": ["this should be delete after without comments"],
              "keep2": ["### this should delete ###", "this should keep"],
              "### delete3 ###": {"this_should_be_delete": "this_should_be_delete"},
              "fix1": {"### delete4 ###": {"this_should_be_delete": []},
                       "keep3":"### this should keep ###",
                       "keep4":["this should be keep",
                                "### this should delete ###",
                                {"### delete5 ###": []}]}},
             "### a fix dict end 1 ###",
             "### another list start 2 ###",
             ["this should be keep"],
             "### another list end 2 ###"]
        self.answer_workflow2_file_path = os.path.join(self.fake_path, "fake_json_with_comments2.json")
        self.answer_workflow2 = \
            {
             "fix_dict": {"name": "first_node",
                          "func": ["None.test.fake.fake_core", "where_am_i", ""],
                          "type": "cold",
                          "inputs": ["input:cmd_params"],
                          "outputs": ["output:results"],
                          "extra_args": [],
                          "extra_kwargs": {},
                          "flags": [],
                          "keep1": "### this should keep ###",
                          "keep2": ["this should keep"],
                          "fix1": {
                                   "keep3": "### this should keep ###",
                                   "keep4": ["this should be keep",
                                             {}]}},
             "a_list": ["this should be keep"],
             }
        self.answer_workflow2_with_comments = \
            {"### a fix dict start ###": "1",
             "fix_dict": {"name": "first_node",
                          "func": ["None.test.fake.fake_core", "where_am_i", ""],
                          "type": "cold",
                          "inputs": ["input:cmd_params"],
                          "outputs": ["output:results"],
                          "extra_args": [],
                          "extra_kwargs": {},
                          "flags": [],
                          "### delete1 ###": "this should be delete after without comments",
                          "keep1": "### this should keep ###",
                          "### delete2 ###": ["this should be delete after without comments"],
                          "keep2": ["### this should delete ###", "this should keep"],
                          "### delete3 ###": {"this_should_be_delete": "this_should_be_delete"},
                          "fix1": {"### delete4 ###": {"this_should_be_delete": []},
                                   "keep3": "### this should keep ###",
                                   "keep4": ["this should be keep",
                                             "### this should delete ###",
                                             {"### delete5 ###": []}]}},
             "### a fix dict end ###": "1",

             "### another list start ###": "2",
             "a_list": ["this should be keep"],
             "### another list end ###": "2"}

    def test_file_to_json(self):
        answer = self.answer_whereami
        json_data = Json.file_to_json(self.whereami_file_path)
        assert json_data == answer

    def test_file_to_json_1(self):
        answer = self.answer_workflow1_with_comments
        json_data = Json.file_to_json(self.answer_workflow1_file_path)
        assert json_data == answer

    def test_file_to_json_without_comments_1(self):
        answer = self.answer_workflow1
        json_data = Json.file_to_json_without_comments(self.answer_workflow1_file_path)
        assert json_data == answer

    def test_file_to_json_2(self):
        answer = self.answer_workflow2_with_comments
        json_data = Json.file_to_json(self.answer_workflow2_file_path)
        assert json_data == answer

    def test_file_to_json_without_comments_2(self):
        answer = self.answer_workflow2
        json_data = Json.file_to_json_without_comments(self.answer_workflow2_file_path)
        assert json_data == answer
