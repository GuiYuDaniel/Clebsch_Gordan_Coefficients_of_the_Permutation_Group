# -*- coding:utf8 -*-
"""
测试cgc_results下所有初始数据库的正确性
"""


import copy
import os
import pytest
import shutil
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from conf.cgc_config import top_path, cgc_rst_folder
from core.cgc_utils.cgc_local_db import CGCLocalDb
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name
from utils.log import get_logger


logger = get_logger(__name__)


class TestDBInit(object):

    def setup_class(self):
        self.s_1 = 1

        self.calculated_tables_s_1 = 1

        self.young_diagrams_s_1 = [[1]]
        _, self.young_diagrams_file_name = get_young_diagrams_file_name(self.s_1)  # S1
        _, self.young_diagrams_full_path = \
            get_young_diagrams_file_name(self.s_1, is_full_path=True)  # <top_path>/cgc_results/young_diagrams_info/S1

    def teardown_class(self):
        pass

    def test_calculated_tables_info(self):
        pass

    def test_young_diagrams_info(self):
        assert os.path.exists(self.young_diagrams_full_path + ".pkl")
        assert os.path.exists(self.young_diagrams_full_path + ".txt")
        flag, data = YoungDiagramInfo(self.s_1).query_by_file_name(self.young_diagrams_file_name)
        assert flag
        assert data
        assert isinstance(data, dict)
        assert data.get("data") == self.young_diagrams_s_1

