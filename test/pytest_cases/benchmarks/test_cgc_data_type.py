# -*- coding:utf8 -*-
"""
测试
更好的储存CGC的数据格式
1，纯字典
2，行列式(sympy.Matrix)（分数）
3，行列式(numpy.array)（整数+N）
"""


import copy
import numpy as np
import sympy as sp
import sys
import os
import pytest
import random
import time
from itertools import product
from conf.cgc_config import cgc_rst_folder, cgc_file_name_format_none
from core.isf_and_cgc import save_cgc
from core.cgc_utils.cgc_local_db import get_cgc_file_name
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


class TestCGCDataType(object):
    """从
    1，大小；
    2，对称性；
    3，读、存时间
    几个方面比较CGC数据格式
    """

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        logger.info("@" * 20)
        logger.info("=" * 20)
        logger.info("[start] TestCGCDataType")

        # 格式div_σ, div_μ, σ_step, μ_step
        self.data_1 = (3, 4, 1, 2)
        self.data_2 = (7, 9, 2, 4)
        self.data_3 = (9, 8, 2, 2)
        self.data_4 = (15, 23, 2, 3)
        self.data_5 = (24, 18, 4, 2)
        self.data_6 = (66, 70, 6, 5)
        self.data_7 = (88, 50, 4, 3)
        self.data_8 = (90, 107, 2, 4)

    def teardown_class(self):
        self.protector.protector_teardown()
        logger.info("[end] TestCGCDataType")
        logger.info("^" * 20)
        logger.info("\n")
        logger.info("\n")

    @staticmethod
    def make_same_cgc_square_fake_data(div_σ, div_μ, σ_step: int, μ_step: int):
        # rst
        fake_cgc_dict = {}  # dict形式的CGC
        fake_cgc_sp_determinant = sp.zeros(div_σ, div_μ)

        # 造fake CGC
        σ_m_list = list(range(1, div_σ + 1, σ_step))
        # μ_m_list = list(range(1, div_μ + 1, μ_step))
        μ_m_list = list(range(div_μ, 0, -μ_step))
        # 分子
        for σ_m, μ_m in product(σ_m_list, μ_m_list):
            random_int = random.randint(1, 20)
            fake_cgc_dict[(σ_m, μ_m)] = random_int
            fake_cgc_sp_determinant[σ_m - 1, μ_m - 1] = random_int
        # 分母
        total = sum(fake_cgc_dict.values())
        for (σ_i, σ_m), (μ_i, μ_m) in product(enumerate(σ_m_list), enumerate(μ_m_list)):
            fake_cgc_dict[(σ_m, μ_m)] = sp.Rational(fake_cgc_dict[(σ_m, μ_m)], total)
        fake_cgc_sp_determinant = fake_cgc_sp_determinant / total

        return fake_cgc_dict, fake_cgc_sp_determinant

    def test_cgc_data_type(self):
        """
        结论：
        从txt上看，越稀疏，字典越节省存储。密度在1/4附近，字典和矩阵占用的存储达到相近水平
        总体而言，pkl化的Matrix还是要稍大一些。
        所以，在未找到更好的办法和字典还能用之前，先继续使用字典。
        """
        for i in range(1, 8 + 1):
            data_tuple = eval("self.data_{}".format(i))
            name_tuple_for_dict = (1, [0], [0], [i], None, 0)
            name_tuple_for_determinant = (1, [0], [0], [i], None, 1)
            _, full_path_for_dict = get_cgc_file_name(*name_tuple_for_dict, is_full_path=True)
            _, full_path_for_determinant = get_cgc_file_name(*name_tuple_for_determinant, is_full_path=True)
            cgc_dict, cgc_sp_determinant = self.make_same_cgc_square_fake_data(*data_tuple)
            flag, msg = save_cgc(*name_tuple_for_dict, cgc_dict)
            assert flag, msg
            flag, msg = save_cgc(*name_tuple_for_determinant, cgc_sp_determinant)
            assert flag, msg

            logger.debug("####")
            logger.debug("for div_σ={}, div_μ={}, σ_step={}, μ_step={}".format(*data_tuple))
            logger.debug("no 0 data number={} in {}, density={}".format(len(cgc_dict), cgc_sp_determinant.shape,
                                                                        len(cgc_dict) / len(cgc_sp_determinant)))
            logger.info("type sp dict used sys.getsizeof={}, os.path.getsize: pkl={}, txt={}".format(
                sys.getsizeof(cgc_dict),
                os.path.getsize(full_path_for_dict + ".pkl"),
                os.path.getsize(full_path_for_dict + ".txt")))
            logger.info("type sp determinant used sys.getsizeof={}, os.path.getsize: pkl={}, txt={}".format(
                sys.getsizeof(cgc_sp_determinant),
                os.path.getsize(full_path_for_determinant + ".pkl"),
                os.path.getsize(full_path_for_determinant + ".txt")))
            logger.debug("\n")
