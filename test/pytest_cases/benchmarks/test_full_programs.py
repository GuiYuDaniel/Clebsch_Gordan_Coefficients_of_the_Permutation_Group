# -*- coding:utf8 -*-
"""
测试
全流程
"""


import os
import time
import math
import numpy as np
import pytest
from conf.cgc_config import cgc_db_name
from core.young_diagrams import create_young_diagrams
from core.branching_laws import create_branching_laws
from core.young_tableaux import create_young_tableaux
from core.yamanouchi_matrix import create_yamanouchi_matrix
from core.characters_and_gi import create_characters_and_gi
from core.cg_series import create_cg_series
from core.eigenvalues import create_eigenvalues
from core.isf_and_cgc import create_isf_and_cgc
from db.local_db_protector import DBProtector
from pytest_cases.benchmarks.utils import get_run_function_n_times as run_func_n_times
from utils.log import get_logger


logger = get_logger(__name__)


# def run_func_n_times(func, *args, n=100, **kwargs):
#     """获取run n次function的时间，内部函数，为了时间更加准确，暂不检查参数了"""
#     start_time = time.time()
#     for i in range(int(n)):
#         func_rst = func(*args, **kwargs)
#     speed_time = time.time() - start_time
#     return "func={} runs {} times with speeding {}s".format(func.__name__, n, speed_time), func_rst


class TestISFAndCGC(object):

    def setup(self):
        self.protector = DBProtector(cgc_db_name, extension_name=".test_isf_and_cgc_protected")
        self.protector.protector_setup()

    def teardown(self):
        self.protector.protector_teardown()
        pass
    
    def _float_rst(self, s_n=2):
        flag, msg = create_young_diagrams(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_branching_laws(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_young_tableaux(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_yamanouchi_matrix(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_characters_and_gi(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_cg_series(s_n)
        assert flag
        assert msg == s_n

        flag, msg = create_eigenvalues(s_n)
        assert flag
        assert msg == s_n

        flag, finish_s_n = create_isf_and_cgc(s_n)
        assert flag
        assert finish_s_n == s_n
        
    def test_pure_float(self):
        """纯float结果运行时间速度"""
        logger.info("==== 对比纯float结果运行时间速度 ====")

        # 1, s_n=3
        timer, _ = run_func_n_times(self._float_rst, n=1, s_n=3)
        logger.warning("纯float timer={} with s_n={}".format(timer, 3))

        # 2, s_n=4
        timer, _ = run_func_n_times(self._float_rst, n=1, s_n=4)
        logger.debug("纯float timer={} with s_n={}".format(timer, 4))
        #
        # 3, s_n=5
        timer, _ = run_func_n_times(self._float_rst, n=1, s_n=5)
        logger.debug("纯float timer={} with s_n={}".format(timer, 5))

        # 4, s_n=6
        timer, _ = run_func_n_times(self._float_rst, n=1, s_n=6)
        logger.debug("纯float timer={} with s_n={}".format(timer, 6))

        # 4, s_n=7
        timer, _ = run_func_n_times(self._float_rst, n=1, s_n=7)
        logger.debug("纯float timer={} with s_n={}".format(timer, 7))
