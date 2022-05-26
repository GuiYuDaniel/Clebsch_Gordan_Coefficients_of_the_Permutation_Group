# -*- coding:utf8 -*-
"""
测试
function concat功能是否正确执行，目录：core/cgc_utils/cgc_utils
"""


from core.cgc_utils.cgc_utils import concat_s_n
from utils.log import get_logger


logger = get_logger(__name__)


class TestConcat(object):

    def test_concat_true(self):
        flag, s_n = concat_s_n(100)
        assert flag
        assert s_n == 100
        flag, s_n = concat_s_n(100, 100)
        assert flag
        assert s_n == 100
        flag, s_n = concat_s_n(a=100)
        assert flag
        assert s_n == 100
        flag, s_n = concat_s_n(a=100, b=100)
        assert flag
        assert s_n == 100
        flag, s_n = concat_s_n(100, a=100)
        assert flag
        assert s_n == 100
        flag, s_n = concat_s_n(100, 100, 100, a=100, b=100)
        assert flag
        assert s_n == 100

    def test_concat_false(self):
        flag, s_n = concat_s_n()
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, 101)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, 101, 102)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, a=101)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(a=100, b=101)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, 101, 100, 100)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, 101, a=100, b=100)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, 100.01)  # 注意：100, 100.0 因为set会把它们合并，是能通过的。但是一般不会出现，就不防这种情况了
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100.0)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
        flag, s_n = concat_s_n(100, s_n=101)
        logger.info("Error is supposed here!")
        assert not flag
        assert isinstance(s_n, str)
