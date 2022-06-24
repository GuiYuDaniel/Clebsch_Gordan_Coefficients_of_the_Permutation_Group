# -*- coding: utf-8 -*-
"""
这段代码用来测试
1，优化前后，isf_matrix_element的一致性，提升性
"""


import copy
import numpy as np
import pytest
import random
import time
from itertools import product
from conf.cgc_config import cgc_rst_folder
from core.cgc_utils.cgc_db_typing import YamanouchiMatrixInfo
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


class S4RealData(object):
    """
    S4取自真实情况  σ μ ν β σ’ μ’ ν’ β’
    # # S4 ISF
    # 取σ=[3, 1] μ=[2, 2]           --> ν可选(3, 1)(2, 1, 1)
    # 取ν=[3, 1] β=1                --> ν’可选(3)(2, 1)
    # 取ν’=[2, 1] β’=1              --> σ’μ’组合有([3][2, 1])和([2, 1][2, 1])
    # 既{((3,), (2, 1), 1): 0.4999999999999999, ((2, 1), (2, 1), 1): 0.4999999999999999}
    # 取σ’μ’=[3][2, 1]
    # 同时，σ’μ’=[2, 1][2, 1]也会被一起算出来
    ##
    # # S3 CGC
    # 取σ’μ’=[3][2, 1] ν’=[2, 1]β’=1 m=1
    # {(1, 1): 1.0, 'N': 1.0}
    # 同时，取σ’μ’=[2, 1][2, 1] ν’=[2, 1]β’=1 m=1
    # {(1, 1): 0.4999999999999999, (2, 2): -0.4999999999999999, 'N': 0.9999999999999998}
    # 转化为S4的partial
    # {(1, 1): 1.0, 'N': 1.0}
    # {(2, 1): 0.4999999999999999, (3, 2): -0.4999999999999999, 'N': 0.9999999999999998}
    ##
    # # S4 IN
    # 取σ=[3, 1]的IN 3*3
    # (14)
    # [[-0.33333333333333315, -0.4714045207910316, -0.8164965809277259],
    #  [-0.4714045207910316, 0.8333333333333331, -0.28867513459481287],
    #  [-0.8164965809277259, -0.28867513459481287, 0.4999999999999999]]
    # (24)
    # [[-0.3333333333333333, -0.4714045207910317, 0.8164965809277259],
    #  [-0.4714045207910317, 0.8333333333333333, 0.28867513459481287],
    #  [0.8164965809277259, 0.28867513459481287, 0.5]]
    # (34)
    # [[-0.3333333333333333, 0.9428090415820634, 0],
    #  [0.9428090415820634, 0.3333333333333333, 0],
    #  [0, 0, 1]]
    # 取μ=[2, 2]的IN 2*2
    # (14)
    # [[-0.4999999999999999, 0.8660254037844386],
    #  [0.8660254037844386, 0.4999999999999999]]
    # (24)
    # [[-0.5, -0.8660254037844386],
    #  [-0.8660254037844386, 0.5]]
    # (34)
    # [[1, 0],
    #  [0, -1]]
    """

    def __init__(self):
        # S4取自真实情况  σ μ ν β σ’ μ’ ν’ β’
        self.s_n, self.matrix_div = [4], [2]
        self.cgc_st_dict_1 = {(1, 1): 1.0, 'N': 1.0}  # 取σ’μ’=[3][2, 1] ν’=[2, 1]β’=1 m=1
        self.cgc_st_dict_2 = {(1, 1): 0.4999999999999999, (2, 2): -0.4999999999999999, 'N': 0.9999999999999998}
        # 将S3的m'，对应到S4的m
        self.cgc_dict_1 = {(1, 1): 1.0, 'N': 1.0}  # 取σ’μ’=[3][2, 1] ν’=[2, 1]β’=1 m=1
        self.cgc_dict_2 = {(2, 1): 0.4999999999999999, (3, 2): -0.4999999999999999, 'N': 0.9999999999999998}  # σ’[2, 1]
        self.in_matrix_1 = {(1, 4): np.array([[-0.33333333333333315, -0.4714045207910316, -0.8164965809277259],
                                              [-0.4714045207910316, 0.8333333333333331, -0.28867513459481287],
                                              [-0.8164965809277259, -0.28867513459481287, 0.4999999999999999]]),
                            (2, 4): np.array([[-0.3333333333333333, -0.4714045207910317, 0.8164965809277259],
                                              [-0.4714045207910317, 0.8333333333333333, 0.28867513459481287],
                                              [0.8164965809277259, 0.28867513459481287, 0.5]]),
                            (3, 4): np.array([[-0.3333333333333333, 0.9428090415820634, 0],
                                              [0.9428090415820634, 0.3333333333333333, 0],
                                              [0, 0, 1]])}
        self.in_matrix_2 = {(1, 4): np.array([[-0.4999999999999999, 0.8660254037844386],
                                              [0.8660254037844386, 0.4999999999999999]]),
                            (2, 4): np.array([[-0.5, -0.8660254037844386],
                                              [-0.8660254037844386, 0.5]]),
                            (3, 4): np.array([[1, 0],
                                              [0, -1]])}
        self.isf_matrix_element_1 = 2
        self.isf_matrix_element_2 = 0


class TestISFMatrixElementOpt(object):
    """优化前后，isf_matrix_element的一致性，提升性

    1，比较计算准确性
    2，比较运算速度
    """

    def setup_class(self):
        self.protector = DBProtector(cgc_rst_folder, extension_name=".test_cgc_db.protector")
        self.protector.protector_setup()

        logger.info("@" * 20)
        logger.info("=" * 20)
        logger.info("[start] TestISFMatrixElementOpt: 优化前后，isf_matrix_element的一致性，提升性")
        # 其他参数
        self.load_db_speed_time = 0.8/10000
        # S4取自真实情况  σ μ ν β σ’ μ’ ν’ β’
        self.s_4_data = S4RealData()
        self.s_n_real, self.matrix_div_real = [4], [2]  # hard code for show
        # 下面的是随机生成的伪数据(s_n不可以有4！！！)
        self.s_n_short, self.matrix_div_short_list = [5, 6], [8, 13]
        self.s_n_long, self.matrix_div_long_list = [7, 8, 9], [27, 50, 100]
        # 用生成时间正比模拟load时间
        self.left_in_matrix_dict, self.right_in_matrix_dict = None, None  # dict(np.ndarray(float))
        self.left_cgc_dict, self.right_cgc_dict = None, None

    def teardown_class(self):
        self.protector.protector_teardown()
        logger.info("[end] TestISFMatrixElementOpt: 优化前后，isf_matrix_element的一致性，提升性")
        logger.info("^" * 20)
        logger.info("\n")
        logger.info("\n")

    def make_fake_data(self, data_flag, s_n, σ, μ, σ_st, μ_st, fake_div=None, fake_density=None):
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, 1, σ, 1, 1,
                                                        fake_div=fake_div, fake_direction="left")
        μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, 1, μ, 1, 1,
                                                        fake_div=fake_div + 3, fake_direction="right")

    def clean_fake_data(self):
        self.left_in_matrix_dict, self.right_in_matrix_dict = None, None  # dict(np.ndarray(float))
        self.left_cgc_dict, self.right_cgc_dict = None, None

    @pytest.mark.skip("pass")
    def test_load_pkl_time(self):
        # save
        fake_finish_s_n = 1000
        fake_nu = [2, 1]
        fake_ij = (2, 3,)
        fake_file_name = "S{}/{}/ij{}".format(fake_finish_s_n, fake_nu, fake_ij)
        fake_table = {
            "file_name": fake_file_name,
            "data": np.array([[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]]),
            "flags": {"speed_time": 0,
                      "total_num": 2}
        }
        fake_table_copy = copy.deepcopy(fake_table)

        db_info = YamanouchiMatrixInfo(fake_finish_s_n)
        flag, msg = db_info.insert(fake_table)
        assert flag
        assert msg is True

        # load
        flag, data = db_info.query_by_file_name(fake_file_name)
        assert flag
        assert np.array_equal(data.get("data"), fake_table_copy.get("data"))

        # benchmark
        # Load func=query_by_file_name runs 10000 times with speeding 0.8292748928070068s
        timer, _ = run_func_n_times(db_info.query_by_file_name, fake_file_name, n=10000)
        logger.debug("# Load {} with data={}".format(timer, "np.array"))

    def _born_cgc_dict(self, data_flag, σ_st, fake_div_σ=None, fake_div_μ=None, fake_direction=None, fake_density=None):
        # 得到的已经是由St m' 升级为Sn m的字典
        # sleep
        time.sleep(self.load_db_speed_time)

        # real data
        if data_flag is True:  # Sn=4
            if σ_st == [3]:  # 本例的区分只有σ_st不同，故只取这一个参数作为判断依据
                return copy.deepcopy(self.s_4_data.cgc_dict_1)
            elif σ_st == [2, 1]:
                return copy.deepcopy(self.s_4_data.cgc_dict_2)
            else:
                raise Exception("input wrong σ_st={}".format(σ_st))

        # fake data
        if fake_direction == "left" and self.left_cgc_dict is not None:
            return copy.deepcopy(self.left_cgc_dict)
        elif fake_direction == "right" and self.right_cgc_dict is not None:
            return copy.deepcopy(self.right_cgc_dict)

        fake_cgc_dict = {}
        σ_m_list = list(range(1, fake_div_σ + 1, fake_density))
        μ_m_list = list(range(1, fake_div_μ + 1, fake_density + 1))
        for σ_m, μ_m in product(σ_m_list, μ_m_list):
            fake_cgc_dict[(σ_m, μ_m)] = random.random()
        fake_cgc_dict["N"] = random.randint(1, fake_div_σ + fake_div_μ)

        if fake_direction == "left":
            self.left_cgc_dict = fake_cgc_dict
            return copy.deepcopy(self.left_cgc_dict)
        elif fake_direction == "right":
            self.right_cgc_dict = fake_cgc_dict
            return copy.deepcopy(self.right_cgc_dict)
        else:
            raise Exception("input wrong fake_direction={}".format(fake_direction))

    def _born_in_matrix_dict(self, data_flag, s_n, i, yd, fake_div=None, fake_direction=None):
        # sleep
        time.sleep(self.load_db_speed_time)

        # real data
        if data_flag is True:  # Sn=4
            if yd == [3, 1]:
                in_matrix_dict = self.s_4_data.in_matrix_1
                return copy.deepcopy(in_matrix_dict.get((i, s_n,)))
            elif yd == [2, 2]:
                in_matrix_dict = self.s_4_data.in_matrix_2
                return copy.deepcopy(in_matrix_dict.get((i, s_n,)))
            else:
                raise Exception("input wrong s_n={}, i={}, yd={}".format(s_n, i, yd))

        # fake data
        if fake_direction == "left" and self.left_in_matrix_dict is not None:
            return copy.deepcopy(self.left_in_matrix_dict.get((i, s_n,)))
        elif fake_direction == "right" and self.right_in_matrix_dict is not None:
            return copy.deepcopy(self.right_in_matrix_dict.get((i, s_n,)))

        fake_in_matrix_dict = {}
        for i in range(1, s_n):
            fake_in_matrix_dict[(i, s_n,)] = np.random.random([fake_div, fake_div])

        if fake_direction == "left":
            self.left_in_matrix_dict = fake_in_matrix_dict
            return copy.deepcopy(self.left_in_matrix_dict.get((i, s_n,)))
        elif fake_direction == "right":
            self.right_in_matrix_dict = fake_in_matrix_dict
            return copy.deepcopy(self.right_in_matrix_dict.get((i, s_n,)))
        else:
            raise Exception("input wrong fake_direction={}".format(fake_direction))

    def _from_in_get_matrix_element(self, data_flag, s_n, i, yd, m_left, m_right, fake_div=None, fake_direction=None):
        start_time = time.time()
        in_matrix = self._born_in_matrix_dict(data_flag, s_n, i, yd, fake_div=fake_div, fake_direction=fake_direction)
        # logger.warning("@@@@ in_matrix={}, m_left={}, m_right={}".format(in_matrix, m_left, m_right))
        in_matrix_element = in_matrix[m_left - 1][m_right - 1]
        # logger.warning("@@@@ single _from_in_get_matrix_element used {}".format(time.time() - start_time))
        return in_matrix_element

    def _calc_isf_matrix_element_origin(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                        fake_div=None, fake_density=None):
        """模拟calc_isf_matrix_element

        Real call 使用必须参数
        P.S. data_flag=True, s_n=4, σ=[3, 1], μ=[2, 2], σ_st=[3] or [2, 1], μ_st=[2, 1]
        Fake call 使用默认参数
        P.S. TODO"""
        matrix_element = 0
        for i in range(1, s_n):
            part_matrix_element = 0
            factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                       fake_div_σ=fake_div,
                                                       fake_div_μ=fake_div + 3,
                                                       fake_direction="left",
                                                       fake_density=fake_density)
            factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                        fake_div_σ=fake_div,
                                                        fake_div_μ=fake_div + 3,
                                                        fake_direction="right",
                                                        fake_density=fake_density + 1)
            left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
            right_n = factor_cgc_right_dict.pop("N")
            for (m_σ_left, m_μ_left), factor_cgc_left in factor_cgc_left_dict.items():
                for (m_σ_right, m_μ_right), factor_cgc_right in factor_cgc_right_dict.items():
                    cgc_left_right = np.sqrt(abs(factor_cgc_left * factor_cgc_right) / (left_n * right_n))
                    cgc_left_right = cgc_left_right * np.sign(factor_cgc_left) * np.sign(factor_cgc_right)
                    σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
                                                                    fake_div=fake_div, fake_direction="left")
                    μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
                                                                    fake_div=fake_div + 3, fake_direction="right")
                    part_matrix_element += cgc_left_right * σ_in_element * μ_in_element
            matrix_element += part_matrix_element
        return matrix_element

    def _calc_isf_matrix_element_opt_with_change_order(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                       fake_div=None, fake_density=None):
        matrix_element = 0
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
        right_n = factor_cgc_right_dict.pop("N")
        for (m_σ_left, m_μ_left), factor_cgc_left in factor_cgc_left_dict.items():
            for (m_σ_right, m_μ_right), factor_cgc_right in factor_cgc_right_dict.items():
                cgc_left_right = np.sqrt(abs(factor_cgc_left * factor_cgc_right) / (left_n * right_n))
                cgc_left_right = cgc_left_right * np.sign(factor_cgc_left) * np.sign(factor_cgc_right)
                in_element_sum = 0
                for i in range(1, s_n):
                    σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
                                                                    fake_div=fake_div, fake_direction="left")
                    μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
                                                                    fake_div=fake_div + 3, fake_direction="right")
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += cgc_left_right * in_element_sum
        return matrix_element

    def _calc_isf_matrix_element_opt_with_change_order_1(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                         fake_div=None, fake_density=None):
        matrix_element = 0
        time_1 = time.time()
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        logger.warning("@@@@ time_1={}".format(time.time() - time_1))
        left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
        right_n = factor_cgc_right_dict.pop("N")
        time_2 = time.time()
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = np.sign(factor_cgc_left) * np.sqrt(abs(factor_cgc_left / left_n))
            left_tmp_dict[left_key] = left_tmp
        right_tmp_dict = {}
        for right_key, factor_cgc_right in factor_cgc_right_dict.items():
            right_tmp = np.sign(factor_cgc_right) * np.sqrt(abs(factor_cgc_right / right_n))
            right_tmp_dict[right_key] = right_tmp
        logger.warning("@@@@ time_2={}".format(time.time() - time_2))
        time_3 = time.time()
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, s_n):
                    σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
                                                                    fake_div=fake_div, fake_direction="left")
                    μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
                                                                    fake_div=fake_div + 3, fake_direction="right")
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += left_tmp * right_tmp * in_element_sum
        logger.warning("@@@@ time_3={}".format(time.time() - time_3))
        return matrix_element

    def _calc_isf_matrix_element_opt_with_change_order_2(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                         fake_div=None, fake_density=None):
        matrix_element = 0
        # time_1 = time.time()
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        # logger.warning("@@@@ time_1={}".format(time.time() - time_1))
        left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
        right_n = factor_cgc_right_dict.pop("N")
        # time_2 = time.time()
        left_tmp_dict = {k: np.sign(v) * np.sqrt(abs(v / left_n)) for k, v in factor_cgc_left_dict.items()}
        right_tmp_dict = {k: np.sign(v) * np.sqrt(abs(v / right_n)) for k, v in factor_cgc_right_dict.items()}
        # logger.warning("@@@@ time_2={}".format(time.time() - time_2))
        # time_3 = time.time()
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, s_n):
                    σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
                                                                    fake_div=fake_div, fake_direction="left")
                    μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
                                                                    fake_div=fake_div + 3, fake_direction="right")
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += left_tmp * right_tmp * in_element_sum
        # logger.warning("@@@@ time_3={}".format(time.time() - time_3))
        return matrix_element

    def _calc_isf_matrix_element_opt_with_change_order_3(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                         fake_div=None, fake_density=None):
        matrix_element = 0
        time_1 = time.time()
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        logger.warning("@@@@ time_1={}".format(time.time() - time_1))
        left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
        right_n = factor_cgc_right_dict.pop("N")
        time_2 = time.time()
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = np.sign(factor_cgc_left) * np.sqrt(abs(factor_cgc_left / left_n))
            left_tmp_dict[left_key] = left_tmp
        right_tmp_dict = {}
        for right_key, factor_cgc_right in factor_cgc_right_dict.items():
            right_tmp = np.sign(factor_cgc_right) * np.sqrt(abs(factor_cgc_right / right_n))
            right_tmp_dict[right_key] = right_tmp
        logger.warning("@@@@ time_2={}".format(time.time() - time_2))
        time_3 = time.time()
        _ = self._from_in_get_matrix_element(data_flag, s_n, 1, σ, 1, 1,
                                             fake_div=fake_div, fake_direction="left")
        _ = self._from_in_get_matrix_element(data_flag, s_n, 1, μ, 1, 1,
                                             fake_div=fake_div + 3, fake_direction="right")
        in_matrix_dict_σ = self.left_in_matrix_dict
        in_matrix_dict_μ = self.right_in_matrix_dict
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, s_n):
                    σ_in_element = in_matrix_dict_σ[(i, s_n)][m_σ_left - 1][m_σ_right - 1]
                    μ_in_element = in_matrix_dict_μ[(i, s_n)][m_μ_left - 1][m_μ_right - 1]
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += left_tmp * right_tmp * in_element_sum
        logger.warning("@@@@ time_3={}".format(time.time() - time_3))
        return matrix_element

    def _calc_isf_matrix_element_opt_with_part_tensor(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                      fake_div=None, fake_density=None):
        time_1 = time.time()
        factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
                                                   fake_div_σ=fake_div,
                                                   fake_div_μ=fake_div + 3,
                                                   fake_direction="left",
                                                   fake_density=fake_density)
        factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
                                                    fake_div_σ=fake_div,
                                                    fake_div_μ=fake_div + 3,
                                                    fake_direction="right",
                                                    fake_density=fake_density + 1)
        logger.warning("@@@@ time_1={}".format(time.time() - time_1))
        left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
        right_n = factor_cgc_right_dict.pop("N")
        time_2 = time.time()
        left_keys_list = list(factor_cgc_left_dict.keys())
        left_values_array = np.array(list(factor_cgc_left_dict.values()))
        right_keys_list = list(factor_cgc_right_dict.keys())
        right_values_array = np.array(list(factor_cgc_right_dict.values()))
        logger.warning("@@@@ time_2={}".format(time.time() - time_2))
        time_3 = time.time()
        left_tmp_array = np.sign(left_values_array) * np.sqrt(np.abs(left_values_array / left_n))
        right_tmp_array = np.sign(right_values_array) * np.sqrt(np.abs(right_values_array / right_n))
        logger.warning("@@@@ time_3={}".format(time.time() - time_3))
        time_4 = time.time()
        left_right_array = np.kron(left_tmp_array, right_tmp_array)
        logger.warning("@@@@ time_4={}".format(time.time() - time_4))

        time_5 = time.time()
        in_matrix_flatten_list = []
        for m_σ_left, m_μ_left in left_keys_list:
            for m_σ_right, m_μ_right in right_keys_list:
                in_element_sum = 0
                for i in range(1, s_n):
                    σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
                                                                    fake_div=fake_div, fake_direction="left")
                    μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
                                                                    fake_div=fake_div + 3, fake_direction="right")
                    in_element_sum += σ_in_element * μ_in_element
                in_matrix_flatten_list.append(in_element_sum)
        logger.warning("@@@@ time_5={}".format(time.time() - time_5))
        in_matrix_flatten_array = np.array(in_matrix_flatten_list)
        time_6 = time.time()
        matrix_element = np.sum(left_right_array * in_matrix_flatten_array)
        logger.warning("@@@@ time_6={}".format(time.time() - time_6))
        return matrix_element

    def _calc_isf_matrix_element_opt_with_matrix_calc(self, data_flag, s_n, σ, μ, σ_st, μ_st,
                                                      fake_div=None, fake_density=None):
        """改完全矩阵乘法好像有点难

        目前四个矩阵，分别是：
                          _        _          _   _
        cgc_left: div: yt(σ') * yt(μ')     -> a * b
        cgc_right: div: yt(σ') * yt(μ')    -> a * b
        SUM(in)σ: div: yt(σ) * yt(σ)       -> A * A
        SUM(in)μ: div: yt(μ) * yt(μ)       -> B * B
        _                 _
        a,a 是A的分支律矩阵; b,b 是B的分支律矩阵
        """
        pass

    def _get_isf_matrix_element_real_answer(self, data_flag, s_n, σ, μ, σ_st, μ_st):
        assert data_flag is True
        assert s_n == 4
        assert σ == [3, 1]
        assert μ == [2, 2]
        assert μ_st == [2, 1]
        if σ_st == [3]:
            return self.s_4_data.isf_matrix_element_1
        elif σ_st == [2, 1]:
            return self.s_4_data.isf_matrix_element_2
        else:
            raise Exception("input wrong s_n={}, σ={}, μ={}, σ_st={}, μ_st={}".format(s_n, σ, μ, σ_st, μ_st))

    @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_real_data(self):
        for s_n, σ_st in product(self.s_n_real, ([3], [2, 1])):
            # data_flag, s_n, σ, μ, σ_st, μ_st
            input_tuple = (True, s_n, [3, 1], [2, 2], σ_st, [2, 1])
            answer = self._get_isf_matrix_element_real_answer(*input_tuple)

            # # benchmark

            # origin
            # origin func=_calc_isf_matrix_element_origin runs 100 times with
            # speeding 0.20842885971069336s with abs(error)=2.220446049250313e-16
            # origin func=_calc_isf_matrix_element_origin runs 100 times with
            # speeding 0.3439159393310547s with abs(error)=0.0
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_origin, *input_tuple, n=100)
            logger.debug("# origin {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

            # change order
            # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 100 times with
            # speeding 0.15840506553649902s with abs(error)=2.220446049250313e-16
            # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 100 times with
            # speeding 0.2922501564025879s with abs(error)=5.551115123125783e-17
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order, *input_tuple, n=100)
            logger.debug("# with_change_order {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

            # change order 1
            # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 100 times with
            # speeding 0.1609480381011963s with abs(error)=2.220446049250313e-16
            # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 100 times with
            # speeding 0.29737401008605957s with abs(error)=0.0
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_1, *input_tuple, n=100)
            logger.debug("# with_change_order_1 {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

            # change order 2
            # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 100 times with
            # speeding 0.1675410270690918s with abs(error)=2.220446049250313e-16
            # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 100 times with
            # speeding 0.2940049171447754s with abs(error)=0.0
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_2, *input_tuple, n=100)
            logger.debug("# with_change_order_2 {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

            # change order 3
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_3, *input_tuple, n=100)
            logger.debug("# with_change_order_3 {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

            # with_part_tensor
            # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 100 times with
            # speeding 0.16979193687438965s with abs(error)=2.220446049250313e-16
            # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 100 times with
            # speeding 0.3040800094604492s with abs(error)=0.0
            timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_part_tensor, *input_tuple, n=100)
            logger.debug("# with_part_tensor {} with abs(error)={}".format(timer, abs(answer - rst)))
            assert np.allclose(rst, answer, atol=0.0000000000001)

    def _isf_matrix_element_opt_by_fake_data(self, s_n, fake_div, fake_density=3, n=10):
        # data_flag, s_n, σ, μ, σ_st, μ_st, fake_div, fake_density
        input_tuple = (False, s_n, None, None, None, None)
        input_dict = {"fake_div": fake_div, "fake_density": fake_density}
        self.make_fake_data(*input_tuple, **input_dict)

        # # benchmark
        # origin
        timer, answer = run_func_n_times(self._calc_isf_matrix_element_origin,
                                         *input_tuple, n=n, **input_dict)
        logger.debug("# origin {} with answer={}".format(timer, answer))

        # # change order
        # timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order,
        #                               *input_tuple, n=n, **input_dict)
        # logger.debug("# with_change_order {} with abs(error)={}".format(timer, abs(answer - rst)))
        # assert np.allclose(rst, answer, atol=0.0000000000001)

        # change order 1
        # timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_1,
        #                               *input_tuple, n=n, **input_dict)
        # logger.debug("# with_change_order_1 {} with abs(error)={}".format(timer, abs(answer - rst)))
        # assert np.allclose(rst, answer, atol=0.0000000000001)

        # change order 2
        # timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_2,
        #                               *input_tuple, n=n, **input_dict)
        # logger.debug("# with_change_order_2 {} with abs(error)={}".format(timer, abs(answer - rst)))
        # assert np.allclose(rst, answer, atol=0.0000000000001)

        # change order 3
        timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_change_order_3,
                                      *input_tuple, n=n, **input_dict)
        logger.debug("# with_change_order_3 {} with abs(error)={}".format(timer, abs(answer - rst)))
        assert np.allclose(rst, answer, atol=0.0000000000001)

        # with_part_tensor
        # timer, rst = run_func_n_times(self._calc_isf_matrix_element_opt_with_part_tensor,
        #                               *input_tuple, n=n, **input_dict)
        # logger.debug("# with_part_tensor {} with abs(error)={}".format(timer, abs(answer - rst)))
        # assert np.allclose(rst, answer, atol=0.0000000000001)

        self.clean_fake_data()

    @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_fake_data_s_n_5_div_8(self):
        # # benchmark rst
        # origin
        # origin func=_calc_isf_matrix_element_origin runs 10 times with speeding
        # 0.5069279670715332s with answer=3.544327168919427
        # change order
        # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 10 times with
        # speeding 0.4842488765716553s with abs(error)=0.0
        # change order 1
        # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 10 timeswith
        # speeding 0.4996480941772461s with abs(error)=0.0
        # change order 2
        # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 10 timeswith
        # speeding 0.502476692199707s with abs(error)=0.0
        # with_part_tensor
        # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 10 times with
        # speeding 0.492156982421875s with abs(error)=0.0
        self._isf_matrix_element_opt_by_fake_data(5, 8, n=10)

    @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_fake_data_s_n_6_div_13(self):
        # # benchmark rst
        # origin
        # origin func=_calc_isf_matrix_element_origin runs 10 times with
        # speeding 3.70307993888855s with answer=22.131690655898574
        # change order
        # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 10 times with
        # speeding 3.6169681549072266s with abs(error)=7.105427357601002e-15
        # change order 1
        # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 10 times with
        # speeding 3.676481008529663s with abs(error)=7.105427357601002e-15
        # change order 2
        # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 10 times with
        # speeding 3.77439284324646s with abs(error)=7.105427357601002e-15
        # with_part_tensor
        # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 10 times with
        # speeding 3.751805305480957s with abs(error)=0.0
        self._isf_matrix_element_opt_by_fake_data(6, 13, n=10)

    @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_fake_data_s_n_7_div_27(self):
        # # benchmark rst
        # origin
        # origin func=_calc_isf_matrix_element_origin runs 5 times with
        # speeding 21.43199586868286s with answer=51.72010653992152
        # change order
        # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 5 times with
        # speeding 21.276501178741455s
        # change order 1
        # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 5 times with
        # speeding 20.873061895370483s with abs(error)=9.947598300641403e-14
        # change order 2
        # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 5 times with
        # speeding 21.16419506072998s with abs(error)=9.947598300641403e-14
        # with_part_tensor
        # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 5 times with
        # speeding 20.983052968978882s with abs(error)=5.684341886080802e-14
        self._isf_matrix_element_opt_by_fake_data(7, 27, n=5)

    @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_fake_data_s_n_8_div_50(self):
        # # benchmark rst
        # origin
        # origin func=_calc_isf_matrix_element_origin runs 1 times with
        # speeding 56.37499690055847s with answer=610.992152049222
        # change order
        # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 1 times with
        # speeding 56.661290884017944s with abs(error)=2.8421709430404007e-12
        # change order 1
        # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 1 times with
        # speeding 55.11927819252014s with abs(error)=2.9558577807620168e-12
        # change order 2
        # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 1 times with
        # speeding 55.25003385543823s with abs(error)=2.9558577807620168e-12
        # change order 3 !!!
        # with_change_order_3 func=_calc_isf_matrix_element_opt_with_change_order_3 runs 1 times with
        # speeding 0.2568497657775879s with abs(error)=7.503331289626658e-12
        # with_part_tensor
        # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 1 times with
        # speeding 55.33107805252075s with abs(error)=2.2737367544323206e-13
        self._isf_matrix_element_opt_by_fake_data(8, 50, n=1)

    # def _calc_isf_matrix_element_opt_with_change_order_1(self, data_flag, s_n, σ, μ, σ_st, μ_st,
    #                                                      fake_div=None, fake_density=None):
    #     matrix_element = 0
    #     time_1 = time.time()
    #     factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
    #                                                fake_div_σ=fake_div,
    #                                                fake_div_μ=fake_div + 3,
    #                                                fake_direction="left",
    #                                                fake_density=fake_density)
    #     factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
    #                                                 fake_div_σ=fake_div,
    #                                                 fake_div_μ=fake_div + 3,
    #                                                 fake_direction="right",
    #                                                 fake_density=fake_density + 1)
    #     logger.warning("@@@@ time_1={}".format(time.time() - time_1))
    #     # time_1=0.0023369789123535156
    #     left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
    #     right_n = factor_cgc_right_dict.pop("N")
    #     time_2 = time.time()
    #     left_tmp_dict = {}
    #     for left_key, factor_cgc_left in factor_cgc_left_dict.items():
    #         left_tmp = np.sign(factor_cgc_left) * np.sqrt(abs(factor_cgc_left / left_n))
    #         left_tmp_dict[left_key] = left_tmp
    #     right_tmp_dict = {}
    #     for right_key, factor_cgc_right in factor_cgc_right_dict.items():
    #         right_tmp = np.sign(factor_cgc_right) * np.sqrt(abs(factor_cgc_right / right_n))
    #         right_tmp_dict[right_key] = right_tmp
    #     logger.warning("@@@@ time_2={}".format(time.time() - time_2))
    #     # time_2=0.001046895980834961
    #     time_3 = time.time()
    #     for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
    #         for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
    #             in_element_sum = 0
    #             for i in range(1, s_n):
    #                 σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
    #                                                                 fake_div=fake_div, fake_direction="left")
    #                 μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
    #                                                                 fake_div=fake_div + 3, fake_direction="right")
    #                 in_element_sum += σ_in_element * μ_in_element
    #             matrix_element += left_tmp * right_tmp * in_element_sum
    #     logger.warning("@@@@ time_3={}".format(time.time() - time_3))
    #     # @@@@ time_3=56.20712184906006
    #     return matrix_element

    # def _calc_isf_matrix_element_opt_with_part_tensor(self, data_flag, s_n, σ, μ, σ_st, μ_st,
    #                                                   fake_div=None, fake_density=None):
    #     time_1 = time.time()
    #     factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
    #                                                fake_div_σ=fake_div,
    #                                                fake_div_μ=fake_div + 3,
    #                                                fake_direction="left",
    #                                                fake_density=fake_density)
    #     factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
    #                                                 fake_div_σ=fake_div,
    #                                                 fake_div_μ=fake_div + 3,
    #                                                 fake_direction="right",
    #                                                 fake_density=fake_density + 1)
    #     logger.warning("@@@@ time_1={}".format(time.time() - time_1))
    #     left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
    #     right_n = factor_cgc_right_dict.pop("N")
    #     time_2 = time.time()
    #     left_keys_list = list(factor_cgc_left_dict.keys())
    #     left_values_array = np.array(list(factor_cgc_left_dict.values()))
    #     right_keys_list = list(factor_cgc_right_dict.keys())
    #     right_values_array = np.array(list(factor_cgc_right_dict.values()))
    #     logger.warning("@@@@ time_2={}".format(time.time() - time_2))
    #     time_3 = time.time()
    #     left_tmp_array = np.sign(left_values_array) * np.sqrt(np.abs(left_values_array / left_n))
    #     right_tmp_array = np.sign(right_values_array) * np.sqrt(np.abs(right_values_array / right_n))
    #     logger.warning("@@@@ time_3={}".format(time.time() - time_3))
    #     time_4 = time.time()
    #     left_right_array = np.kron(left_tmp_array, right_tmp_array)
    #     logger.warning("@@@@ time_4={}".format(time.time() - time_4))
    #
    #     time_5 = time.time()
    #     in_matrix_flatten_list = []
    #     for m_σ_left, m_μ_left in left_keys_list:
    #         for m_σ_right, m_μ_right in right_keys_list:
    #             in_element_sum = 0
    #             for i in range(1, s_n):
    #                 σ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, σ, m_σ_left, m_σ_right,
    #                                                                 fake_div=fake_div, fake_direction="left")
    #                 μ_in_element = self._from_in_get_matrix_element(data_flag, s_n, i, μ, m_μ_left, m_μ_right,
    #                                                                 fake_div=fake_div + 3, fake_direction="right")
    #                 in_element_sum += σ_in_element * μ_in_element
    #             in_matrix_flatten_list.append(in_element_sum)
    #     logger.warning("@@@@ time_5={}".format(time.time() - time_5))
    #     in_matrix_flatten_array = np.array(in_matrix_flatten_list)
    #     time_6 = time.time()
    #     matrix_element = np.sum(left_right_array * in_matrix_flatten_array)
    #     logger.warning("@@@@ time_6={}".format(time.time() - time_6))
    #     return matrix_element

    # def _calc_isf_matrix_element_opt_with_change_order_3(self, data_flag, s_n, σ, μ, σ_st, μ_st,
    #                                                      fake_div=None, fake_density=None):
    #     matrix_element = 0
    #     time_1 = time.time()
    #     # time_1=0.0022788047790527344
    #     factor_cgc_left_dict = self._born_cgc_dict(data_flag, σ_st,
    #                                                fake_div_σ=fake_div,
    #                                                fake_div_μ=fake_div + 3,
    #                                                fake_direction="left",
    #                                                fake_density=fake_density)
    #     factor_cgc_right_dict = self._born_cgc_dict(data_flag, μ_st,
    #                                                 fake_div_σ=fake_div,
    #                                                 fake_div_μ=fake_div + 3,
    #                                                 fake_direction="right",
    #                                                 fake_density=fake_density + 1)
    #     logger.warning("@@@@ time_1={}".format(time.time() - time_1))
    #     left_n = factor_cgc_left_dict.pop("N")  # pop的返回值就是删掉的value
    #     right_n = factor_cgc_right_dict.pop("N")
    #     time_2 = time.time()
    #     left_tmp_dict = {}
    #     for left_key, factor_cgc_left in factor_cgc_left_dict.items():
    #         left_tmp = np.sign(factor_cgc_left) * np.sqrt(abs(factor_cgc_left / left_n))
    #         left_tmp_dict[left_key] = left_tmp
    #     right_tmp_dict = {}
    #     for right_key, factor_cgc_right in factor_cgc_right_dict.items():
    #         right_tmp = np.sign(factor_cgc_right) * np.sqrt(abs(factor_cgc_right / right_n))
    #         right_tmp_dict[right_key] = right_tmp
    #     logger.warning("@@@@ time_2={}".format(time.time() - time_2))
    #     # time_2=0.0010459423065185547
    #     time_3 = time.time()
    #     _ = self._from_in_get_matrix_element(data_flag, s_n, 1, σ, 1, 1,
    #                                          fake_div=fake_div, fake_direction="left")
    #     _ = self._from_in_get_matrix_element(data_flag, s_n, 1, μ, 1, 1,
    #                                          fake_div=fake_div + 3, fake_direction="right")
    #     in_matrix_dict_σ = self.left_in_matrix_dict
    #     in_matrix_dict_μ = self.right_in_matrix_dict
    #     for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
    #         for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
    #             in_element_sum = 0
    #             for i in range(1, s_n):
    #                 σ_in_element = in_matrix_dict_σ[(i, s_n)][m_σ_left - 1][m_σ_right - 1]
    #                 μ_in_element = in_matrix_dict_μ[(i, s_n)][m_μ_left - 1][m_μ_right - 1]
    #                 in_element_sum += σ_in_element * μ_in_element
    #             matrix_element += left_tmp * right_tmp * in_element_sum
    #     logger.warning("@@@@ time_3={}".format(time.time() - time_3))
    #     # time_3=0.2528800964355469
    #     return matrix_element

    # @pytest.mark.skip("pass")
    def test_isf_matrix_element_opt_by_fake_data_s_n_9_div_100(self):
        """结论：order 3可用"""
        # # benchmark rst
        # origin
        # origin func=_calc_isf_matrix_element_origin runs 1 times with
        # speeding 901.2666981220245s with answer=2872.689945573261
        # change order
        # with_change_order func=_calc_isf_matrix_element_opt_with_change_order runs 1 times with
        # speeding 901.7362079620361s with abs(error)=1.0186340659856796e-10
        # change order 1
        # with_change_order_1 func=_calc_isf_matrix_element_opt_with_change_order_1 runs 1 times with
        # speeding 899.528153181076s with abs(error)=1.0186340659856796e-10
        # change order 2
        # with_change_order_2 func=_calc_isf_matrix_element_opt_with_change_order_2 runs 1 times with
        # speeding 899.4670078754425s with abs(error)=1.0186340659856796e-10
        # change order 3 !!!
        # with_change_order_3 func=_calc_isf_matrix_element_opt_with_change_order_3 runs 1 times with
        # speeding 3.832380771636963s with abs(error)=1.1823431123048067e-11
        # with_part_tensor
        # with_part_tensor func=_calc_isf_matrix_element_opt_with_part_tensor runs 1 times with
        # speeding 889.4094867706299s with abs(error)=3.637978807091713e-12
        self._isf_matrix_element_opt_by_fake_data(9, 100, n=1)
