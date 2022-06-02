# -*- coding: utf-8 -*-
"""
这段代码用来测试
1，纯numpy PK numpy+fractions.Fractions PK sympy
2，Fraction PK Rational
"""


import time
import numpy as np
import sympy
import pytest
# from sympy.core.numbers import Half
from sympy import Rational, Matrix
from fractions import Fraction
from utils.log import get_logger


logger = get_logger(__name__)


def _get_run_function_n_times(func, *args, n=100, **kwargs):
    """获取run n次function的时间，内部函数，为了时间更加准确，暂不检查参数了"""
    start_time = time.time()
    for i in range(int(n)):
        func_rst = func(*args, **kwargs)
    speed_time = time.time() - start_time
    return "func={} runs {} times with speeding {}s".format(func.__name__, n, speed_time), func_rst


class TestOne(object):
    """纯numpy PK numpy+fractions.Fractions PK sympy

    1，比较装载速度
    2，比较运算速度：
        a）矩阵求和
        b）矩阵直和
        c）矩阵点乘
        d）矩阵本征值
    3，TODO 比较save/load速度
    """

    def setup_class(self):
        logger.info("@" * 20)
        logger.info("=" * 20)
        logger.info("[start] TestOne: 纯numpy PK numpy+fractions.Fractions PK sympy")
        self.data_100_100_list = [[i + j for j in range(100)] for i in range(100)]
        self.data_500_500_list = [[i + j for j in range(500)] for i in range(500)]

    def teardown_class(self):
        logger.info("[end] TestOne: 纯numpy PK numpy+fractions.Fractions PK sympy")
        logger.info("^" * 20)
        logger.info("\n")
        logger.info("\n")

    @pytest.mark.skip("pass")
    def test_make_matrix(self):
        """对比装载速度"""
        logger.info("==== 对比装载速度 ====")

        # 1, 纯numpy
        # func=array runs 100 times with speeding 0.06806182861328125s
        timer, _ = _get_run_function_n_times(np.array, self.data_100_100_list)
        logger.debug("纯numpy {} with data={}".format(timer, "self.data_100_100_list"))

        # func=array runs 100 times with speeding 1.5680510997772217
        timer, _ = _get_run_function_n_times(np.array, self.data_500_500_list)
        logger.debug("纯numpy {} with data={}".format(timer, "self.data_500_500_list"))

        # 2，纯sympy
        # MutableDenseMatrix runs 100 times with speeding 1.8198168277740479s
        timer, _ = _get_run_function_n_times(sympy.Matrix, self.data_100_100_list)
        logger.debug("纯sympy {} with data={}".format(timer, "self.data_100_100_list"))

        # MutableDenseMatrix runs 100 times with speeding 44.70412993431091s
        # 500差不多是100的 5**2 倍  到这里已经很难使用了
        timer, _ = _get_run_function_n_times(sympy.Matrix, self.data_500_500_list)
        logger.debug("纯sympy {} with data={}".format(timer, "self.data_500_500_list"))

    def test_matrix_four_fundamental(self):
        """对比四则运算速度"""
        pass


class TestTwo(object):
    """Fraction PK Rational

    1，比较单纯类的装载速度
    2，比较单纯类的计算速度：
        a）加减
        b）乘除
    3，TODO 比较save/load速度
    """

    def setup_class(self):
        logger.info("@" * 20)
        logger.info("=" * 20)
        logger.info("[start] TestTwo: Fraction PK Rational")
        self.data_100_100_list = [[i + j for j in range(100)] for i in range(100)]

    def teardown_class(self):
        logger.info("[end] TestTwo: Fraction PK Rational")
        logger.info("^" * 20)
        logger.info("\n")
        logger.info("\n")

    @staticmethod
    def _make_method_list(method, size):
        if method == float:
            return [[i / (j + 1) for j in range(size)] for i in range(size)]
        else:
            return [[method(i, j + 1) for j in range(size)] for i in range(size)]

    def test_001_mix_fraction_and_rational(self):
        a = Fraction(1, 3)
        b = Rational(1, 3)
        assert isinstance(a, Fraction)
        assert isinstance(b, Rational)
        aa = np.array([[a for j in range(3)] for i in range(3)])
        bb = np.array([[b for j in range(3)] for i in range(3)])

        # 1， 简单四则，会得到Rational
        c = a + b
        d = b + a
        assert isinstance(a, Fraction)
        assert isinstance(b, Rational)
        assert isinstance(c, Rational) and isinstance(d, Rational)
        c = a * b
        d = b * a
        assert isinstance(c, Rational) and isinstance(d, Rational)

        # 2， np加 乘
        cc = aa + bb  # same as np.add
        dd = bb + aa  # same as np.add
        assert isinstance(cc[0][0], Rational) and isinstance(dd[0][0], Rational)
        cc = aa * bb
        dd = bb * aa
        assert isinstance(cc[0][0], Rational) and isinstance(dd[0][0], Rational)
        cc = aa @ bb
        dd = bb @ aa
        assert isinstance(cc[0][0], Rational) and isinstance(dd[0][0], Rational)

    @pytest.mark.skip("pass")
    def test_make_data(self):
        """对比装载速度"""
        '''
        结论：单纯列表解析，用Fraction是 完全可以接受的 ；至于Rational，需要看真实时间尺度
             ^^^^^^^^^    ^^^^^^^^  ^^^^^^^^^^^
        1， 装载单个分数，Fraction比Rational慢4倍左右，但速度都在接受范围内
        2， 用在列表解析上，Fraction反超Rational4倍多，幸好性能几乎线性（且对比自己，仍然符合维度指数倍，线性于数据量）
        （将method去掉，只构建i, j + 1元组，也就比Fraction快10倍，这个开销还是可以接受的）
        3， 由列表装载到array上，几乎无差别
        '''
        logger.info("==== 对比装载速度 ====")

        # 1, Fraction
        # Fraction func=Fraction runs 100000 times with speeding 0.12648391723632812s
        timer, _ = _get_run_function_n_times(Fraction, 1, 10, n=100000)
        logger.debug("Fraction {} with data={}".format(timer, "1/10"))

        # Fraction List func=_make_method_list runs 100 times with speeding 1.6068532466888428s
        timer, fraction_list_100 = _get_run_function_n_times(self._make_method_list, Fraction, 100)
        logger.debug("Fraction List {} with data={}".format(timer, "100_100"))
        # Fraction List func=_make_method_list runs 10 times with speeding 4.391110897064209s
        timer, fraction_list_500 = _get_run_function_n_times(self._make_method_list, Fraction, 500, n=10)
        logger.debug("Fraction List {} with data={}".format(timer, "500_500"))

        # Fraction Array func=array runs 100 times with speeding 1.078965187072754s
        timer, _ = _get_run_function_n_times(np.array, fraction_list_100)
        logger.debug("Fraction Array {} with data={}".format(timer, "100_100"))
        # Fraction Array func=array runs 10 times with speeding 2.702277898788452s
        timer, _ = _get_run_function_n_times(np.array, fraction_list_500, n=10)
        logger.debug("Fraction Array {} with data={}".format(timer, "500_500"))

        # 2, Rational
        # Rational func=Rational runs 100000 times with speeding 0.030122995376586914s
        timer, _ = _get_run_function_n_times(Rational, 1, 10, n=100000)
        logger.debug("Rational {} with data={}".format(timer, "1/10"))

        # Rational List func=_make_method_list runs 100 times with speeding 7.454913139343262s
        timer, rational_list_100 = _get_run_function_n_times(self._make_method_list, Rational, 100)
        logger.debug("Rational List {} with data={}".format(timer, "100_100"))
        # Rational List func=_make_method_list runs 10 times with speeding 19.09254813194275s
        timer, rational_list_500 = _get_run_function_n_times(self._make_method_list, Rational, 500, n=10)
        logger.debug("Rational List {} with data={}".format(timer, "500_500"))

        # Rational Array func=array runs 100 times with speeding 1.071699857711792s
        timer, _ = _get_run_function_n_times(np.array, rational_list_100)
        logger.debug("Rational Array {} with data={}".format(timer, "100_100"))
        # Rational Array func=array runs 10 times with speeding 2.685981273651123s
        timer, _ = _get_run_function_n_times(np.array, rational_list_500, n=10)
        logger.debug("Rational Array {} with data={}".format(timer, "500_500"))

        # assert
        for fraction_row, rational_row in zip(fraction_list_100, rational_list_100):
            for i, j in zip(fraction_row, rational_row):
                assert str(i) == str(j)

    @pytest.mark.skip("pass")
    def test_matrix_four_fundamental(self):
        """对比四则运算速度"""
        '''
        结论：这回的时间就比较尴尬了，Rational+sympy，忘了它吧；Fraction，加法还行，尽量别乘
        
        1， 相加，Fraction比Rational快3倍多，但仍然比float慢三个量级
        2， 点乘，Fraction比Rational快3倍多，已经比float慢3到4个量级了
        3， sympy自带的Matrix还不如numpy的Array呢
        '''

        float_array_100 = np.array(self._make_method_list(float, 100))
        float_array_200 = np.array(self._make_method_list(float, 200))
        fraction_array_100 = np.array(self._make_method_list(Fraction, 100))
        fraction_array_200 = np.array(self._make_method_list(Fraction, 200))
        rational_array_100 = np.array(self._make_method_list(Rational, 100))
        rational_array_200 = np.array(self._make_method_list(Rational, 200))
        rational_matrix_100 = Matrix(self._make_method_list(Rational, 100))
        rational_matrix_200 = Matrix(self._make_method_list(Rational, 200))

        # 0, Base float
        # Float Add func=add runs 20 times with speeding 0.00019621849060058594s
        timer, _ = _get_run_function_n_times(np.add, float_array_100, float_array_100, n=20)
        logger.debug("Float Add {} with data={}".format(timer, "100_100"))
        # Float Add func=add runs 2 times with speeding 0.0002841949462890625s
        timer, _ = _get_run_function_n_times(np.add, float_array_200, float_array_200, n=2)
        logger.debug("Float Add {} with data={}".format(timer, "200_200"))

        # Float Dot func=dot runs 10 times with speeding 0.004316091537475586s
        timer, _ = _get_run_function_n_times(np.dot, float_array_100, float_array_100, n=10)
        logger.debug("Float Dot {} with data={}".format(timer, "100_100"))
        # Float Dot func=dot runs 10 times with speeding 0.0023148059844970703s
        timer, _ = _get_run_function_n_times(np.dot, float_array_200, float_array_200, n=10)
        logger.debug("Float Dot {} with data={}".format(timer, "200_200"))

        # 1，Fraction
        # Fraction Add func=add runs 20 times with speeding 0.5159449577331543s
        timer, _ = _get_run_function_n_times(np.add, fraction_array_100, fraction_array_100, n=20)
        logger.debug("Fraction Add {} with data={}".format(timer, "100_100"))
        # Fraction Add func=add runs 2 times with speeding 0.2199268341064453s
        timer, _ = _get_run_function_n_times(np.add, fraction_array_200, fraction_array_200, n=2)
        logger.debug("Fraction Add {} with data={}".format(timer, "200_200"))

        # Fraction Dot func=dot runs 1 times with speeding 5.43730092048645s
        timer, _ = _get_run_function_n_times(np.dot, fraction_array_100, fraction_array_100, n=1)
        logger.debug("Fraction Dot {} with data={}".format(timer, "100_100"))
        # Fraction Dot func=dot runs 1 times with speeding 50.333991289138794s
        timer, _ = _get_run_function_n_times(np.dot, fraction_array_200, fraction_array_200, n=1)
        logger.debug("Fraction Dot {} with data={}".format(timer, "200_200"))

        # 2，Rational Array
        # Rational Add func=add runs 20 times with speeding 1.8065309524536133s
        timer, _ = _get_run_function_n_times(np.add, rational_array_100, rational_array_100, n=20)
        logger.debug("Rational Add {} with data={}".format(timer, "100_100"))
        # Rational Add func=add runs 2 times with speeding 0.8531489372253418s
        timer, _ = _get_run_function_n_times(np.add, rational_array_200, rational_array_200, n=2)
        logger.debug("Rational Add {} with data={}".format(timer, "200_200"))

        # Rational Dot func=dot runs 1 times with speeding 21.420315742492676s
        timer, _ = _get_run_function_n_times(np.dot, rational_array_100, rational_array_100, n=1)
        logger.debug("Rational Dot {} with data={}".format(timer, "100_100"))
        # Rational Dot func=dot runs 1 times with speeding 195.694757938385s
        timer, _ = _get_run_function_n_times(np.dot, rational_array_200, rational_array_200, n=1)
        logger.debug("Rational Dot {} with data={}".format(timer, "200_200"))

        # 3，Rational Matrix
        def _matrix_add(a, b):
            return a + b

        def _matrix_dot(a, b):
            return a * b

        # Rational Matrix Add func=_matrix_add runs 20 times with speeding 1.9100239276885986s
        timer, _ = _get_run_function_n_times(_matrix_add, rational_matrix_100, rational_matrix_100, n=20)
        logger.debug("# Rational Matrix Add {} with data={}".format(timer, "100_100"))
        # Rational Matrix Add func=_matrix_add runs 2 times with speeding 0.8602540493011475s
        timer, _ = _get_run_function_n_times(_matrix_add, rational_matrix_200, rational_matrix_200, n=2)
        logger.debug("# Rational Matrix Add {} with data={}".format(timer, "200_200"))

        # Rational Matrix Dot func=_matrix_dot runs 1 times with speeding 23.167457818984985s
        timer, _ = _get_run_function_n_times(_matrix_dot, rational_matrix_100, rational_matrix_100, n=1)
        logger.debug("# Rational Matrix Dot {} with data={}".format(timer, "100_100"))
        # Rational Matrix Dot func=_matrix_dot runs 1 times with speeding 215.38406085968018s
        timer, _ = _get_run_function_n_times(_matrix_dot, rational_matrix_200, rational_matrix_200, n=1)
        logger.debug("# Rational Matrix Dot {} with data={}".format(timer, "200_200"))
