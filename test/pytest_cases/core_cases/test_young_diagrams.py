# -*- coding:utf8 -*-
"""
测试
core/young_diagrams.py
下所有功能是否正确执行
"""


import pytest
import time
from conf.cgc_config import default_s_n
from core.young_diagrams import calc_single_young_diagrams
from utils.log import get_logger


logger = get_logger(__name__)


class TestYoungDiagrams(object):

    def setup_class(self):
        self.young_diagrams_s_1 = [[1]]
        self.young_diagrams_s_2 = [[2], [1, 1]]
        self.young_diagrams_s_3 = [[3], [2, 1], [1, 1, 1]]
        self.young_diagrams_s_4 = [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        self.young_diagrams_s_5 = [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        self.young_diagrams_s_6 = [[6], [5, 1],
                                   [4, 2], [4, 1, 1],
                                   [3, 3], [3, 2, 1], [3, 1, 1, 1],
                                   [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1, 1]]
        self.young_diagrams_s_7 = [[7], [6, 1],
                                   [5, 2], [5, 1, 1],
                                   [4, 3], [4, 2, 1], [4, 1, 1, 1],
                                   [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [3, 1, 1, 1, 1],
                                   [2, 2, 2, 1], [2, 2, 1, 1, 1], [2, 1, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1, 1, 1]]
        self.default_s_n = default_s_n

    def teardown_class(self):
        pass

    def test_calc_single_young_diagrams_sample(self):
        head_s_n = 1
        tail_s_n = 7
        head_time = time.time()
        for s_n in range(head_s_n, tail_s_n + 1):
            flag, young_diagrams = calc_single_young_diagrams(s_n, recursion_deep=7)
            logger.debug("for s_n={}, test young_diagrams={}".format(s_n, young_diagrams))
            assert flag
            assert young_diagrams == eval("self.young_diagrams_s_{}".format(s_n)), "s_n={}".format(s_n)
        tail_time = time.time()
        logger.debug("calc young_diagrams from S{} to S{}, used {}s".format(head_s_n, tail_s_n, tail_time - head_time))

    def test_calc_single_young_diagrams_out_of_recursion_deep(self):
        s_n = 1000
        flag, young_diagrams = calc_single_young_diagrams(s_n, recursion_deep=7)
        logger.debug("for s_n={}, test young_diagrams={}".format(s_n, young_diagrams))
        assert flag
        assert young_diagrams is False

    @pytest.mark.skip("pass")
    def test_calc_single_young_diagrams_recursion(self):
        pass
