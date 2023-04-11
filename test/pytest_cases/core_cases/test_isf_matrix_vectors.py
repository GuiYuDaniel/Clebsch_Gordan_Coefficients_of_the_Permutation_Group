# -*- coding:utf8 -*-


import time
import sympy as sp
from sympy import Rational as Ra
from sympy import sqrt, Matrix
from itertools import chain
import pytest
from core.young_diagrams import load_young_diagrams
from core.young_tableaux import load_young_table_num
from core.eigenvalues import is_eigenvalue, is_eigenvector
from core.symmetry_combination import calc_ϵ_map_dicts
from core.isf_and_cgc import ΣMDataHelper, DataHelper, CalcHelper, ISFHelper, CGCHelper, load_ϵ, EHelper
from utils.log import get_logger


logger = get_logger(__name__)


class Data(object):

    def __init__(self):
        # 测试参数
        # self.st = 6
        # self.sn = 7
        # self.σ = [5, 1, 1]
        # self.μ = [4, 2, 1]
        # self.ν_st = [3, 1, 1, 1]
        #
        # self.st = 2
        # self.sn = 3
        # self.σ = [2, 1]
        # self.μ = [2, 1]
        # self.ν_st = [2]
        # sq = [(Ra(1)/2, -Ra(1)/2)]
        # # self.ν_st = [1, 1]
        # # sq = [(-Ra(1)/2, -Ra(1)/2)]
        # self.ν = [2, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 3
        # self.sn = 4
        # self.σ = [3, 1]
        # self.μ = [3, 1]
        # self.ν_st = [2, 1]
        # sq = [(1,)]
        # self.ν = [2, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 3
        # self.sn = 4
        # self.σ = [3, 1]
        # self.μ = [2, 1, 1]
        # self.ν_st = [3]
        # sq = [(1,)]
        # self.ν = [3, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 4
        # self.sn = 5
        # self.σ = [3, 1, 1]
        # self.μ = [3, 1, 1]
        # self.ν_st = [3, 1]
        # sq = [(0, Ra(1)/2, -Ra(1)/2, 0)]
        # # sq = [(sp.symbols("sq1 sq2 sq3 sq4"))]
        # # self.ν_st = [2, 1, 1]
        # # sq = [(Ra(1)/2, 0, 0, Ra(1)/2)]
        # self.ν = [3, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        # '''
        # 让[311]*[311]=[311]对称不一致的：
        # νμσ, σνμ, νσμ, μνσ, ν~μ~σ, σ~ν~μ, ν~σ~μ, μ~ν~σ, ν~μσ~, σ~νμ~, ν~σμ~, μ~νσ~, νμ~σ~, σν~μ~, νσ~μ~, μν~σ~
        # '''

        # self.st = 4
        # self.sn = 5
        # self.σ = [3, 1, 1]
        # self.μ = [3, 1, 1]
        # self.ν_st = [3, 1]
        # sq = [(Ra(1)/2, 0, 0, -Ra(1)/2),
        #       (Ra(1)/12, Ra(5)/12, Ra(5)/12, Ra(1)/12)]
        # # self.ν_st = [2, 1, 1]
        # # sq = [(Ra(1)/2, 0, 0, Ra(1)/2)]
        # self.ν = [3, 2]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 4
        # self.sn = 5
        # self.σ = [3, 1, 1]
        # self.μ = [3, 1, 1]
        # self.ν_st = [2, 2]
        # sq = [(Ra(5)/16, Ra(3)/16, Ra(3)/16, -Ra(5)/16),
        #       (0, Ra(1)/2, -Ra(1)/2, 0)]
        # # self.ν_st = [2, 1, 1]
        # # sq = [(0, -Ra(1)/2, -Ra(1)/2, 0),
        # #       (-Ra(5)/12, Ra(1)/12, -Ra(1)/12, Ra(5)/12)]
        # self.ν = [2, 2, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [4, 2]
        # self.μ = [4, 1, 1]
        # # self.ν_st = [4, 1]
        # # sq = [(Ra(1)/2, Ra(1)/12, -Ra(1)/4, Ra(1)/6),
        # #       (Ra(5)/18, -Ra(5)/12, Ra(49)/180, Ra(1)/30)]
        # self.ν_st = [3, 1, 1]
        # sq = [(-Ra(1)/18, -Ra(8)/27, Ra(1)/9, Ra(1)/2, -Ra(1)/27),
        #       (Ra(5)/18, 0, Ra(1)/45, Ra(1)/10, Ra(3)/5)]
        # self.ν = [4, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [4, 2]
        # self.μ = [3, 2, 1]
        # # self.ν_st = [3, 2]
        # # sq = [(0, Ra(16)/75, 0, Ra(18)/25, -Ra(1)/15, 0),
        # #       (Ra(4)/9, -Ra(2)/15, 0, Ra(1)/180, -Ra(1)/6, Ra(1)/4),
        # #       (Ra(5)/72, -Ra(1)/300, -Ra(5)/8, Ra(8)/225, Ra(4)/15, 0)]
        # # self.ν_st = [3, 1, 1]
        # # sq = [(Ra(8)/45, Ra(64)/135, Ra(8)/45, Ra(1)/18, 0, Ra(8)/135, Ra(1)/18),
        # #       (-Ra(1)/9, 0, Ra(1)/9, Ra(5)/36, Ra(1)/2, 0, -Ra(5)/36),
        # #       (-Ra(1)/360, Ra(1)/60, -Ra(1)/360, -Ra(2)/9, 0, Ra(8)/15, -Ra(2)/9)]
        # self.ν_st = [2, 2, 1]
        # sq = [(0, Ra(16)/75, 0, 0, Ra(1)/15, Ra(18)/25),
        #       (0, Ra(2)/15, -Ra(4)/9, Ra(1)/4, -Ra(1)/6, -Ra(1)/180),
        #       (Ra(5)/8, -Ra(1)/300, Ra(5)/72, 0, -Ra(4)/15, Ra(8)/225)]
        # self.ν = [3, 2, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [4, 2]
        # self.μ = [3, 2, 1]
        # self.ν_st = [3, 1, 1]
        # sq = [(-Ra(25)/72, Ra(25)/108, -Ra(1)/72, Ra(5)/72, -Ra(1)/4, -Ra(1)/54, -Ra(5)/72),
        #       (Ra(5)/144, -Ra(5)/216, -Ra(5)/16, Ra(4)/9, 0, Ra(5)/27, 0)]
        # self.ν = [3, 1, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [4, 1, 1]
        # self.μ = [4, 1, 1]
        # # self.ν_st = [4, 1]
        # # sq = [(Ra(9)/20, -Ra(3)/40, -Ra(3)/40, -Ra(2)/5),
        # #       (Ra(1)/4, Ra(3)/8, Ra(3)/8, 0)]
        # self.ν_st = [3, 2]
        # sq = [(-Ra(9)/50, Ra(3)/25, Ra(3)/25, Ra(27)/50, -Ra(1)/25),
        #       (Ra(49)/250, Ra(3)/125, Ra(3)/125, Ra(27)/250, Ra(81)/125)]
        # self.ν = [4, 2]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [4, 1, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [3, 2]
        # sq = [(Ra(3)/8, Ra(1)/20, Ra(3)/200, 0, -Ra(2)/5, 0, Ra(4)/25),
        #       (0, Ra(3)/400, -Ra(5)/8, Ra(3)/20, -Ra(3)/200, Ra(81)/400, 0),
        #       (0, Ra(8)/25, 0, -Ra(2)/5, Ra(1)/25, Ra(6)/25, 0),
        #       (0, Ra(5)/16, Ra(3)/200, Ra(1)/4, Ra(9)/40, -Ra(3)/80, Ra(4)/25)]
        # # self.ν_st = [3, 1, 1]
        # # sq = [(Ra(1)/24, -Ra(1)/4, Ra(1)/24, -Ra(1)/3, 0, 0, -Ra(1)/3, 0),
        # #       (-Ra(1)/160, 0, Ra(1)/160, Ra(1)/80, -Ra(27)/160, -Ra(5)/8, -Ra(1)/80, Ra(27)/160),
        # #       (-Ra(4)/15, 0, -Ra(4)/15, -Ra(1)/30, -Ra(1)/5, 0,  -Ra(1)/30, -Ra(1)/5),
        # #       (-Ra(25)/96, 0, Ra(25)/96, -Ra(3)/16, Ra(1)/32, -Ra(1)/24, Ra(3)/16, -Ra(1)/32)]
        # # self.ν_st = [2, 2, 1]
        # # sq = [(-Ra(3)/200, Ra(1)/20, Ra(3)/8, Ra(4)/25, -Ra(2)/5, 0, 0),
        # #       (-Ra(5)/8, -Ra(3)/400, 0, 0, Ra(3)/200, Ra(81)/400, Ra(3)/20),
        # #       (0, Ra(8)/25, 0, 0, Ra(1)/25, -Ra(6)/25, Ra(2)/5),
        # #       (Ra(3)/200, -Ra(5)/16, 0, -Ra(4)/25, -Ra(9)/40, -Ra(3)/80, Ra(1)/4)]
        # self.ν = [3, 2, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 3]
        # self.μ = [3, 2, 1]
        # self.ν_st = [3, 2]
        # sq = [(Ra(3)/4, Ra(1)/10, -Ra(3)/20),  # 书中解
        #       (0, Ra(3)/5, Ra(2)/5)]
        # # sq = [(Ra(3)/28, Ra(7)/10, Ra(27)/140),  # 我的解
        # #       (Ra(9)/14, 0, -Ra(5)/14)]
        # # self.ν_st = [3, 1, 1]
        # # sq = [(-Ra(1)/12, Ra(5)/6, 0, Ra(1)/12),
        # #       (Ra(1)/2, 0, 0, Ra(1)/2)]
        # # self.ν_st = [2, 2, 1]
        # # # sq = [(-Ra(3)/20, Ra(1)/10, Ra(3)/4),  # 书中解
        # # #       (-Ra(2)/5, -Ra(3)/5, 0)]
        # # sq = [(-Ra(15)/28, -Ra(5)/14, -Ra(3)/28),  # 我的解
        # #       (-Ra(1)/70, Ra(12)/35, -Ra(9)/14)]
        # self.ν = [3, 2, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        self.st = 5
        self.sn = 6
        self.σ = [3, 2, 1]
        self.μ = [3, 2, 1]
        self.ν_st = [3, 2]
        sq = [(Ra(3)/40, -Ra(1)/16, 0, -Ra(1)/16, 0, -Ra(3)/5, -Ra(1)/16, 0, -Ra(1)/16, Ra(3)/40),
              (0, -Ra(1)/4, 0, Ra(1)/4, 0, 0, -Ra(1)/4, 0, Ra(1)/4, 0),
              (0, -Ra(1)/16, -Ra(3)/8, Ra(1)/16, 0, 0, Ra(1)/16, Ra(3)/8, -Ra(1)/16, 0),
              (Ra(9)/40, 0, Ra(1)/8, 0, -Ra(3)/10, 0, 0, Ra(1)/8, 0, -Ra(9)/40),
              (Ra(3)/80, Ra(1)/8, -Ra(3)/16, Ra(1)/8, -Ra(1)/20, 0, -Ra(1)/8, -Ra(3)/16, -Ra(1)/8, -Ra(3)/80)]
        # sq = [(Ra(9)/40, 0, Ra(1)/8, 0, -Ra(3)/10, 0, 0, Ra(1)/8, 0, -Ra(9)/40),
        #       (0, -Ra(1)/4, 0, Ra(1)/4, 0, 0, -Ra(1)/4, 0, Ra(1)/4, 0)]
        # sq = [(0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
        #       (0, 0, 0, 0, 2, 0, 0, 0, 0, 0),
        #       (0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
        #       (0, 0, 0, 0, 0, 2, 0, 0, 0, 0),
        #       (0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
        #       (0, 0, 0, 0, 2, 2, 0, 0, 0, 0)]
        # # sq = [(Ra(3)/40, -Ra(1)/16, 0, -Ra(1)/16, 0, -Ra(3)/5, -Ra(1)/16, 0, -Ra(1)/16, Ra(3)/40)]
        # sq = [(Ra(9)/40, 0, Ra(1)/8, 0, -Ra(3)/10, 0, 0, Ra(1)/8, 0, -Ra(9)/40),
        #       (Ra(3)/80, Ra(1)/8, -Ra(3)/16, Ra(1)/8, -Ra(1)/20, 0, -Ra(1)/8, -Ra(3)/16, -Ra(1)/8, -Ra(3)/80)]
        # # sq_pro = [(sp.symbols("τ1c0 τ1c1 τ1c2 τ1c3 τ1c4 τ1c5 τ1c6 τ1c7 τ1c8 τ1c9")),
        # #           (sp.symbols("τ2c0 τ2c1 τ2c2 τ2c3 τ2c4 τ2c5 τ2c6 τ2c7 τ2c8 τ2c9"))]
        # s_soe = [sq[0], sq[3], sq[4]]
        # a_soe = [sq[1], sq[2]]
        # # xa, xb, xc, xd, xe = sp.symbols("xa xb xc xd xe")
        # # # # sq = [(sp.sign(s_soe[0][i]) * sp.sqrt(abs(s_soe[0][i])) * xa +
        # # # #        sp.sign(s_soe[1][i]) * sp.sqrt(abs(s_soe[1][i])) * xb +
        # # # #        sp.sign(s_soe[2][i]) * sp.sqrt(abs(s_soe[2][i])) * xc for i in range(10))]
        # # # # sq = [(sp.sign(s_soe[0][i]) * sp.sqrt(abs(s_soe[0][i])) * -2 * xc +
        # # # #        sp.sign(s_soe[1][i]) * sp.sqrt(abs(s_soe[1][i])) * xb +
        # # # #        sp.sign(s_soe[2][i]) * sp.sqrt(abs(s_soe[2][i])) * xc for i in range(10))]
        # sq = [(sp.sign(s_soe[0][i]) * sp.sqrt(abs(s_soe[0][i])) * (sqrt(2)) +
        #        sp.sign(s_soe[2][i]) * sp.sqrt(abs(s_soe[2][i])) * 1 for i in range(10))]
        # sq = [(-Ra(1)/80, Ra(1)/6, -Ra(1)/16, Ra(1)/6, -Ra(1)/60, Ra(2)/5, 0, -Ra(1)/16, 0, -Ra(9)/80),
        #       (Ra(9) / 40, 0, Ra(1) / 8, 0, -Ra(3) / 10, 0, 0, Ra(1) / 8, 0, -Ra(9) / 40)]
        # self.ν_st = [3, 1, 1]
        # sq = [(0, -Ra(1)/48, -Ra(1)/8, -Ra(5)/24, -Ra(1)/48, -Ra(1)/8,
        #        0, -Ra(1)/48, -Ra(1)/8, -Ra(5)/24, -Ra(1)/48, -Ra(1)/8, 0),
        #       (-Ra(5)/24, 0, 0, Ra(5)/24, 0, 0, -Ra(1)/6, 0, 0, -Ra(5)/24, 0, 0, -Ra(5)/24),
        #       (Ra(5)/24, -Ra(1)/48, Ra(1)/8, 0, Ra(1)/48, -Ra(1)/8,
        #        0, Ra(1)/48, -Ra(1)/8, 0, -Ra(1)/48, Ra(1)/8, -Ra(5)/24),
        #       (0, -Ra(1)/4, 0, 0, -Ra(1)/4, 0, 0, Ra(1)/4, 0, 0, Ra(1)/4, 0, 0),
        #       (0, 0, -Ra(1)/4, 0, 0, -Ra(1)/4, 0, 0, Ra(1)/4, 0, 0, Ra(1)/4, 0)]
        # self.ν_st = [2, 2, 1]
        self.ν = [3, 2, 1]
        self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        # self.sym_key = "σ~μ~ν"

        # 书中没有的
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [5]
        # sq = []
        # self.ν_st = [4, 1]
        # self.ν = [5, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [4, 1]
        # sq = []
        # # self.ν_st = [3, 2]
        # self.ν = [4, 2]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [4, 1]
        # sq = []
        # # self.ν_st = [3, 1, 1]
        # self.ν = [4, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [3, 2]
        # sq = []
        # self.ν = [3, 3]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [3, 1, 1]
        # sq = []
        # self.ν_st = [2, 1, 1, 1]
        # self.ν = [3, 1, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [2, 2, 1]
        # sq = []
        # self.ν = [2, 2, 2]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [2, 2, 1]
        # sq = []
        # self.ν_st = [2, 1, 1, 1]
        # self.ν = [2, 2, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν
        #
        # self.st = 5
        # self.sn = 6
        # self.σ = [3, 2, 1]
        # self.μ = [3, 2, 1]
        # self.ν_st = [2, 1, 1, 1]
        # sq = []
        # self.ν_st = [1, 1, 1, 1, 1]
        # self.ν = [2, 1, 1, 1, 1]
        # self.sym_key = "νμσ"  # 一个指定的对称组合，可以得到一组确定的sym_σμν

        # 防bug
        '''
        这个bug来自sympy1.0版本之后的类注册器。如果不手动建立一次非整数的Rational，让它的注册表里存在Rational，直接读pickle会报错。

        C, including its class ClassRegistry, has been deprecated since SymPy
        1.0. It will be last supported in SymPy version 1.0. Use direct
        imports from the defining module instead. See
        https://github.com/sympy/sympy/issues/9371 for more info.
        '''
        _ = Matrix([[Ra(1) / 2]])
        # 需要加工的数据
        self.vector_list = [Matrix([sp.sign(j) * sqrt(abs(j)) for j in i]) for i in sq]
        # self.tmp_vector_list = [Matrix([j for j in i]) for i in sq]
        # self.vector_list = sp.GramSchmidt(self.tmp_vector_list, True)


class TestISFMatrix(object):

    def setup_class(self):
        # 它不改变数据库，所以不需要保护
        # 需要提前准备好S6的全部数据库

        self.data = Data()
        self.st = self.data.st
        self.sn = self.data.sn
        self.σ = self.data.σ
        self.μ = self.data.μ
        self.ν_st = self.data.ν_st
        self.ν = self.data.ν
        self.sym_key = self.data.sym_key
        self.vector_list = self.data.vector_list
        logger.info("Sn={}: σ={}, μ={}, ν_st={}, ν={}".format(self.sn, self.σ, self.μ, self.ν_st, self.ν))

        # 准备
        ϵ_key2groups_dict, groups2ϵ_key_dict = calc_ϵ_map_dicts()

        self.data_st = DataHelper(self.st, ϵ_key2groups_dict, groups2ϵ_key_dict)
        data_st_st_yt_num_list = []
        _, yd_st_st_list = load_young_diagrams(self.sn - 1 - 1)
        for yd_st_st in yd_st_st_list:
            _, total_num = load_young_table_num(self.sn - 1 - 1, yd_st_st)
            data_st_st_yt_num_list.append(total_num)
        self.data_st.add_more_data(data_st_st_yt_num_list, yd_st_st_list)

        self.data_sn = DataHelper(self.sn, ϵ_key2groups_dict, groups2ϵ_key_dict)
        self.data_sn.add_more_data(self.data_st.yt_num_list, self.data_st.yd_list)

        self.data_σ_μ = ΣMDataHelper(self.sn, self.σ, self.μ, self.data_sn, self.data_st)
        self.data_σ_μ.hook_meta_ν_list_of_σμ_and_get_ν_st_list_inner_meta_bl()

        self.isf_func = ISFHelper()
        self.isf_func.enable_now_s_n(self.sn)

        self.ϵ_func = EHelper()
        self.ϵ_func.enable_now_s_n(self.sn)

    def teardown_class(self):
        pass

    @pytest.mark.skip("pass")
    def test_isf_matrix(self):
        """测试多种计算isf_matrix方法，是否能得出一致的结果"""
        # 准备
        row_index_list = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        assert row_index_list, \
            "row_index_list={} must real list, pls check with data_σ_μ={}, ν_st={}, data_st={}" \
            "".format(row_index_list, self.data_σ_μ, self.ν_st, self.data_st)
        logger.info("row_index_list={}".format(row_index_list))

        # 测试
        # 4_193a
        start_time_4_193a = time.time()
        flag, isf_matrix_with_4_193a = \
            self.isf_func._calc_isf_matrix_4_193a(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        assert flag, \
            "_calc_isf_matrix_4_193a by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail with msg={}" \
            "".format(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, isf_matrix_with_4_193a)
        speed_time_4_193a = int(time.time() - start_time_4_193a)
        logger.info("isf_matrix_with_4_193a={}".format(isf_matrix_with_4_193a))

        # last_m
        start_time_last_m = time.time()
        flag, isf_matrix_with_last_m = \
            self.isf_func._calc_isf_matrix_with_last_m(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        assert flag, \
            "_calc_isf_matrix_with_last_m by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail with msg={}" \
            "".format(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, isf_matrix_with_4_193a)
        speed_time_last_m = int(time.time() - start_time_last_m)
        logger.info("isf_matrix_with_last_m={}".format(isf_matrix_with_last_m))

        # first_m
        start_time_first_m = time.time()
        flag, isf_matrix_with_first_m = \
            self.isf_func._calc_isf_matrix_with_first_m(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ)
        assert flag, \
            "_calc_isf_matrix_with_first_m by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail with msg={}" \
            "".format(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, isf_matrix_with_4_193a)
        speed_time_first_m = int(time.time() - start_time_first_m)
        logger.info("isf_matrix_with_first_m={}".format(isf_matrix_with_first_m))

        # time
        if self.sn >= 6 and speed_time_4_193a != 0 and speed_time_last_m != 0 and speed_time_first_m != 0:
            logger.debug("isf_matrix_with_4_193a using time={}\n"
                         "isf_matrix_with_last_m using time={}\n"
                         "isf_matrix_with_first_m using time={}\n"
                         "".format(speed_time_4_193a, speed_time_last_m, speed_time_first_m))

    @pytest.mark.skip("pass")
    def test_isf_eigenvector(self):
        """测试书中数据，是否是isf_matrix的本征矢量"""
        # 准备
        row_index_list = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        flag, isf_matrix = \
            self.isf_func._calc_isf_matrix_4_193a(row_index_list, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        logger.info("isf_matrix={}".format(isf_matrix))

        # 测试
        rst_list = []
        for single_i, single_vector in enumerate(self.vector_list):
            flag, is_ev = is_eigenvector(single_vector, isf_matrix)
            assert flag, "err_msg={}".format(is_ev)
            rst_list.append(is_ev)
            if is_ev is False:
                logger.error("single_i={}, single_vector={}".format(single_i, single_vector))
        assert all(rst_list), "{}".format(rst_list)

        '''
        小结1：(τ默认等于2)
        a1，[4, 2] * [4, 1, 1] = [4, 1, 1](by[4, 1])先开始，书中的ISF也是本征矢
        a2，[4, 2] * [4, 1, 1] = [4, 1, 1](by[3, 1, 1])先开始，书中的ISF也是本征矢
        
        b1，[4, 2] * [3, 2, 1] = [3, 2, 1](by[3, 2])先开始，书中的ISF也是本征矢(τ=3)
        b2，[4, 2] * [3, 2, 1] = [3, 2, 1](by[3, 1, 1])先开始，书中的ISF也是本征矢(τ=3)
        b3，[4, 2] * [3, 2, 1] = [3, 2, 1](by[2, 2, 1])先开始，书中的ISF也是本征矢(τ=3)
        
        c1，[4, 1, 1] * [4, 1, 1] = [4, 2](by[4, 1])先开始，书中的ISF也是本征矢
        c2，[4, 1, 1] * [4, 1, 1] = [4, 2](by[3, 2])先开始，书中的ISF也是本征矢
        
        d1，[4, 1, 1] * [3, 2, 1] = [3, 2, 1](by[3, 2])先开始，书中的ISF也是本征矢(τ=4)
        d2，[4, 1, 1] * [3, 2, 1] = [3, 2, 1](by[3, 1, 1])先开始，书中的ISF 不 不 不 是本征矢!(τ=4)（第一个是，其他不是）
        d3，[4, 1, 1] * [3, 2, 1] = [3, 2, 1](by[2, 2, 1])先开始，书中的ISF也是本征矢(τ=4)
        
        e1，[3, 3] * [3, 2, 1] = [3, 2, 1](by[3, 2])先开始，书中的ISF也是本征矢
        e2，[3, 3] * [3, 2, 1] = [3, 2, 1](by[3, 1, 1])先开始，书中的ISF也是本征矢
        e3，[3, 3] * [3, 2, 1] = [3, 2, 1](by[2, 2, 1])先开始，书中的ISF 不 不 不 是本征矢!（第一个不是，第二个是）
        
        f1，[3, 2, 1] * [3, 2, 1] = [3, 2, 1](by[3, 2])先开始，书中的ISF也是本征矢(τ=5)
        f2，[3, 2, 1] * [3, 2, 1] = [3, 2, 1](by[3, 1, 1])先开始，书中的ISF 不 不 不 是本征矢!(τ=5)（第4是，其他不是）
        
        小结2：
        几种猜想：
        1，4-195c不总能得到满足eigenvector的ISF
        2，单纯书错了
        3，我的也可以
        4，其实，存在既满足4-195c，也满足ev的特殊线性组合，前人没有意识到，需要搞定
        '''

    # @pytest.mark.skip("pass")
    def test_isf_4_195c_eigenvector(self):
        """
        测试由公式4-195c（英文书4-195b）和绝对相位ISF计算出来的相对ISF，是否是相应isf_matrix的本征矢量
        1，是否一定是isf_matrix的eigenvector（结论：是！）
        2，书中答案是否相互自洽（结论：不是！）
        """
        # 准备
        # 公共部分
        ν_bl_yds = self.data_sn.get_bl_yds(self.ν)  # ν_st的list
        m_ν_st_st = 1
        # 绝对相位的部分
        ab_ν_st = ν_bl_yds[0]
        ab_row_index_list = self.isf_func.calc_row_indexes(self.data_σ_μ, ab_ν_st, self.data_st)
        ab_vector_list = self.vector_list
        ab_isf_square_list = [Matrix([sp.sign(j) * (j ** 2) for j in i]) for i in self.vector_list]
        logger.info("ab_ν_st={}".format(ab_ν_st))
        logger.info("ab_row_index_list={}".format(ab_row_index_list))
        logger.info("ab_isf_square_list={}".format(ab_isf_square_list))
        # # 测试本征性
        flag, isf_matrix = \
            self.isf_func._calc_isf_matrix_4_193a(ab_row_index_list, ab_ν_st, self.data_sn, self.data_σ_μ,
                                                  self.data_st)
        # rst_list = []
        # for single_i, single_vector in enumerate(ab_vector_list):
        #     flag, is_ev = is_eigenvector(single_vector, isf_matrix)
        #     assert flag, "err_msg={}".format(is_ev)
        #     rst_list.append(is_ev)
        #     if is_ev is False:
        #         logger.error("single_i={}, single_vector={}".format(single_i, single_vector))
        # assert all(rst_list), "{}".format(rst_list)

        # 相对相位的部分
        re_ν_st_list = ν_bl_yds[1::]
        for re_ν_st in re_ν_st_list:
            # 下面三个m的意义分别是：
            # 如此做的根本原因是，去掉Sn和Sn-1，ab_ν和re_ν它们的ν_st_st一致，且m''一致
            # 1 m_ν：               令当前ν_st（非ν的第一分支）的m'取1时，ν的m；
            # 2 m_ν_st_fbl：        令当前ν_st的第一分支ν_st_st的m''取1时，ν_st_fbl（ν的第一分支）的m'
            #                       (ν_st_st是ν_st的第一分支，可以证明它也是ν_st_fbl的分支，但不一定是第一分支)
            # 3,m_ν_by_m_ν_st_fbl： 按照2中方法，确定m_ν_st_fbl后（也意味着确定了ν_st_fbl的杨盘），
            #                       ν对应的m
            ν_st_st_of_re_ν_st = self.data_st.get_bl_yds(re_ν_st)[0]  # 它也必然是ab_ν_st的分支
            m_st__re_ν_st_of_ν = self.data_st.quick_calc_m(m_ν_st_st, ν_st_st_of_re_ν_st, re_ν_st)  # 必是1
            assert m_st__re_ν_st_of_ν == 1, "m_st__re_ν_st_of_ν={} should be 1 but not".format(m_st__re_ν_st_of_ν)
            m__re_ν = self.data_sn.quick_calc_m(m_st__re_ν_st_of_ν, re_ν_st, self.ν)  # 1 m_ν
            m_st__ab_ν_st_of_re = self.data_st.quick_calc_m(m_ν_st_st, ν_st_st_of_re_ν_st, ab_ν_st)  # 2 m_ν_st_fbl
            m__ab_ν_of_ab_ν_st_of_re = \
                self.data_sn.quick_calc_m(m_st__ab_ν_st_of_re, ab_ν_st, self.ν)  # 3,m_ν_by_m_ν_st_fbl

            re_row_index_list = self.isf_func.calc_row_indexes(self.data_σ_μ, re_ν_st, self.data_st)

            # 测试
            # 测试由绝对相位ISF生成的相对相位ISF的计算
            phase_vector_list = []
            phase_square_list = []
            # 使用4-195c正算phase_vector，放弃seo_vector
            for ab_isf_vector in ab_vector_list:  # 按τ展开
                phase_vector_element_list = []
                for re_row in re_row_index_list:
                    flag, single_phase_vector_new = \
                        self.isf_func._calc_isf_by_known_isf(ab_isf_vector, ab_row_index_list, ab_ν_st,
                                                             m_st__ab_ν_st_of_re, m__ab_ν_of_ab_ν_st_of_re,
                                                             re_row, re_ν_st, m_st__re_ν_st_of_ν, m__re_ν,
                                                             self.ν, self.data_sn, self.data_σ_μ)
                    assert flag
                    phase_vector_element_list.append(single_phase_vector_new)
                phase_vector = Matrix(phase_vector_element_list)
                phase_vector_list.append(phase_vector)
                phase_square = Matrix([sp.sign(i) * (i**2) for i in phase_vector_element_list])
                phase_square_list.append(phase_square)
            logger.info("re_ν_st={}".format(re_ν_st))
            logger.info("re_row_index_list={}".format(re_row_index_list))
            logger.info("re_isf_square_list={}".format(phase_square_list))
            # # 测试本征性
            # flag, isf_matrix = \
            #     self.isf_func._calc_isf_matrix_4_193a(re_row_index_list, re_ν_st, self.data_sn, self.data_σ_μ,
            #                                           self.data_st)
            # rst_list = []
            # for single_i, single_vector in enumerate(phase_vector_list):
            #     flag, is_ev = is_eigenvector(single_vector, isf_matrix)
            #     assert flag, "err_msg={}".format(is_ev)
            #     rst_list.append(is_ev)
            #     if is_ev is False:
            #         logger.error("single_i={}, single_vector={}".format(single_i, single_vector))
            # assert all(rst_list), "{}".format(rst_list)

        '''
        小结1：(τ默认等于2)
        a, [4, 2] * [4, 1, 1] = [4, 1, 1]
        [4, 1](首)分支：抄书上的，是本征矢
        [3, 1, 1]分支：是本征矢，且和书一致
        
        b, [4, 2] * [3, 2, 1] = [3, 2, 1](τ=3)
        [3, 2](首)分支：抄书上的，是本征矢
        [3, 1, 1]分支：是本征矢，且和书一致
        [2, 2, 1]分支：是本征矢，且和书一致
        
        c, [4, 1, 1] * [4, 1, 1] = [4, 2]
        [4, 1](首)分支：抄书上的，是本征矢
        [3, 2]分支：是本征矢，且和书一致
        
        d, [4, 1, 1] * [3, 2, 1] = [3, 2, 1](τ=4)
        [3, 2](首)分支：抄书上的，是本征矢
        [3, 1, 1]分支：绝对值能对上，但符号有的对不上（完全一致的只有第1个，因为是0）(与书相比，最后一行全部反号！)；但是4-195c出来的却是本征矢
        上面这个错误的"特点"是，它们符号相反的行标都为([3, 1, 1], [2, 2, 1], 2)
        [2, 2, 1]分支：是本征矢，且和书一致
        
        e, [3, 3] * [3, 2, 1] = [3, 2, 1]
        [3, 2](首)分支：抄书上的，是本征矢
        [3, 1, 1]分支：是本征矢，且和书一致
        [2, 2, 1]分支：是本征矢，且和书一致
        
        f，[3, 2, 1] * [3, 2, 1] = [3, 2, 1](τ=5)
        [3, 2](首)分支：抄书上的，是本征矢
        [3, 1, 1]分支：绝对值能对上，但符号有的对不上（完全一致的只有第4个，因为是0）(与书相比，最后一行全部反号！)；但是4-195c出来的却是本征矢
        上面这个错误的"特点"是，它们符号相反的行标都为([3, 1, 1], [3, 1, 1])和([2, 2, 1], [3, 1, 1], 2)
        
        g, [3, 3] * [3, 2, 1] = [3, 1, 1, 1]
        [3, 1, 1](首)分支：抄书上的，是本征矢
        [2, 1, 1, 1]分支：是本征矢，且和书一致
        
        
        小结2：
        对比test_isf_eigenvector和test_isf_sym_eigenvector可以看出：
        1，4-195c是满足本征性的
        2，或者书错了，或者我的isf_matrix（Sn-1的结论）错了
        3，书中数据，至少比我的差，因为它不满足本征性了
        
        
        小结3：
        也就是说，我只要能完美得到一个分支的ISF，就万事大吉了！
        '''

    @pytest.mark.skip("pass")
    def test_isf_calc_soe_and_sym_vector(self):
        """
        测试按照下面两种顺序计算的isf和对称isf是否正确，也检查抛弃的soe是否和书中一致
        1，先计算A*A=B构型的ISF，再根据对称性计算A*B=A构型的ISF，检查它们是否都正确
        2，先计算A*B=A构型的ISF，再根据对称性计算A*A=B构型的ISF，检查它们是否都正确
        """
        # 指定的对称
        meta_σμν = (self.σ, self.μ, self.ν)
        sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(self.sym_key)
        sym_σμν = tuple(meta_σμν[d] if k is False else
                        self.data_sn.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
        sym_σ, sym_μ, sym_ν = sym_σμν
        # sym_ν_st = [4, 1]
        logger.debug("meta:{}*{}={} by {}\n"
                     "sym:{}*{}={}"
                     "".format(self.σ, self.μ, self.ν, self.ν_st, sym_σ, sym_μ, sym_ν))
        # 计算meta_rows和sym_rows（对称的rows可能对应不同isf）
        meta_rows = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        sym_list = []  # (sym_row, sym_ν_st)
        for single_meta_row in meta_rows:
            meta_row_σ, meta_row_μ, meta_row_τ = single_meta_row if len(single_meta_row) == 3 \
                else (*single_meta_row, None)
            meta_σμν_st = (meta_row_σ, meta_row_μ, self.ν_st)
            sym_σμν_st = tuple(meta_σμν_st[d] if k is False else
                               self.data_st.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
            single_sym_row = (sym_σμν_st[0], sym_σμν_st[1]) if meta_row_τ is None \
                else (sym_σμν_st[0], sym_σμν_st[1], meta_row_τ)
            sym_list_element = (single_sym_row, sym_σμν_st[2])
            sym_list.append(sym_list_element)
        logger.debug("meta_rows={}, sym_list={}".format(meta_rows, sym_list))

        flag, isf_matrix = \
            self.isf_func._calc_isf_matrix_4_193a(meta_rows, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        logger.info("isf_matrix={}".format(isf_matrix))

        # 修改_calc_schmidt_orthogonalization_tricky源码，得到指定算法的soe，并测试
        λ_ν_st = self.data_st.get_eigenvalue(self.ν_st)
        λ_ν = self.data_sn.get_eigenvalue(self.ν)
        e_value = λ_ν - λ_ν_st
        e_vectors = self.isf_func._calc_eigenspace_by_eigenvalue(isf_matrix, e_value)
        τ_max = len(e_vectors)

        logger.debug("e_vectors={}, τ_max={}".format(e_vectors, τ_max))

        soe_vectors = self.isf_func._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if τ_max > 1 \
            else [self.isf_func._calc_orthogonalization_vector(e_vectors[0])]

        logger.debug("soe_vectors={}".format(soe_vectors))

        # 这里尝试用其他的候选soe计算

        sym_σμν = (sym_σ, sym_μ, sym_ν)
        sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(self.sym_key)

        helper = Helper(self.st, self.sn)

        for i, this_τ_soe_vector in enumerate(soe_vectors):
            sym_τ = i + 1
            this_τ_meta_isf_square = [sp.sign(i) * (i**2) for i in this_τ_soe_vector]
            logger.info("\n{}*{}={}τ{} by {} ==> {}*{}={}τ{}, rows type below:"
                        "this_τ_meta_isf_square={}"
                        "".format(self.σ, self.μ, self.ν, sym_τ, self.ν_st, sym_σ, sym_μ, sym_ν, sym_τ,
                                  this_τ_meta_isf_square))

            logger_msg = "{:<5}, {:<10} {:<30} = {:<15} from {:<30} = {:<16}\n" \
                         "".format("Index", "sym_ν_st", "sym_row", "sym_isf_element",
                                   "meta_row", "meta_isf_element")
            for (meta_i, meta_row), s_element in zip(enumerate(meta_rows), sym_list):
                sym_row, sym_ν_st = s_element
                sym_isf_element = helper.calc_sym_isf_element(sym_row, sym_σμν, sym_τ, sym_ν_st, self.sym_key,
                                                              sym_d3, sym_k4,
                                                              self.data_sn, self.data_st,
                                                              this_τ_meta_isf_square[meta_i])
                logger_msg += "{:<5}, {:<10} {:<30} = {:<15} from {:<30} = {:<16}\n" \
                              "".format(meta_i, str(sym_ν_st), str(sym_row), str(sym_isf_element),
                                        str(meta_row), str(this_τ_meta_isf_square[meta_i]))
            logger.info("\n{}\n".format(logger_msg))

        '''
        小结1：
        a, 从 [4, 2] * [4, 1, 1] = [4, 1, 1] 到 [4, 1, 1] * [4, 1, 1] = [4, 2]
        [4, 1](首)分支：本身解与书中不一致，本身解的对称与书中也不一致；
                      抛弃解与书中不一致，抛弃解的对称与书中也不一致
        [3, 1, 1]分支：本身解与书中不一致，本身解的对称与书中也不一致；
                      抛弃解与书中不一致，抛弃解的对称与书中也不一致
                      
        b, 从 [4, 2] * [3, 2, 1] = [3, 2, 1](τ=3) 到 [3, 2, 1] * [3, 2, 1] = [4, 2](τ=3)
        [3, 2](首)分支：本身解只有1/3组一致，对称解无答案；对称解row的对称性为：sxx
        [3, 1, 1]分支：pass
        [2, 2, 1]分支：pass
        
        c, 从 [4, 1, 1] * [4, 1, 1] = [4, 2] 到 [4, 2] * [4, 1, 1] = [4, 1, 1]
        # [4, 1](首)分支：本身解与书中一致!本身解的对称与书中也一致!!!
        [3, 2]分支：本身解与书中不一致，本身解的对称与书中也不一致；
                   抛弃解与书中不一致，抛弃解的对称与书中也不一致
                   
        d, 从 [4, 1, 1] * [3, 2, 1] = [3, 2, 1](τ=4) 到 [3, 2, 1] * [3, 2, 1] = [4, 1, 1](τ=4)
        [3, 2](首)分支：啥也不是
        [3, 1, 1]分支：pass
        [2, 2, 1]分支：pass
        
        e, 从 [3, 3] * [3, 2, 1] = [3, 2, 1] 到 [3, 2, 1] * [3, 2, 1] = [3, 3]
        [3, 2](首)分支：啥也不是
        [3, 1, 1]分支：本身解与书中一致!本身解的对称与书中也一致!!!
        [2, 2, 1]分支：本身解与书中一致!本身解的对称与书中也一致!!!
        
        
        B, 从 [3, 2, 1] * [3, 2, 1] = [4, 2](τ=3) 到 [4, 2] * [3, 2, 1] = [3, 2, 1](τ=3)
        # [4, 1](首)分支：本身解与书中一致!本身解的对称与书中也一致!!!
        [3, 2]分支：本身解只有1/3组一致
        [3, 1, 1]分支：pass
        
        C, 从 [3, 2, 1] * [3, 2, 1] = [4, 1, 1](τ=4) 到 [4, 1, 1] * [3, 2, 1] = [3, 2, 1](τ=4)
        [4, 1](首)分支：本身解只有1/4组一致，对称性是saaa
        [3, 1, 1]分支：本身解只有1/4组一致，对称性是(s+a)xxx
        
        D, 从 [3, 2, 1] * [3, 2, 1] = [3, 3] 到 [3, 3] * [3, 2, 1] = [3, 2, 1]
        # [3, 2](首)分支：本身解与书中一致!本身解的对称与书中也一致!!!
        
        
        小结2：
        看起来，如果先做AAA，再做AAB，能得到一致解。其他解，用对称和4-195c是可以完成的。
        '''

    @pytest.mark.skip("pass")
    def test_isf_calc_main_vector(self):
        """
        这里，都是AAA和AAB形式的内积。测试目标：
        1，检查τ=2的，首分支是否都是正确soe
        2，τ>2的，想办法让首分支与书中的线性组合相同
        """
        # 指定的对称
        # meta_σμν = (self.σ, self.μ, self.ν)
        # sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(self.sym_key)
        # sym_σμν = tuple(meta_σμν[d] if k is False else
        #                 self.data_sn.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
        # sym_σ, sym_μ, sym_ν = sym_σμν
        # # sym_ν_st = [4, 1]
        # logger.debug("meta:{}*{}={} by {}\n"
        #              "sym:{}*{}={}"
        #              "".format(self.σ, self.μ, self.ν, self.ν_st, sym_σ, sym_μ, sym_ν))
        # 计算meta_rows和sym_rows（对称的rows可能对应不同isf）
        meta_rows = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        flag, isf_matrix = \
            self.isf_func._calc_isf_matrix_4_193a(meta_rows, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        logger.info("isf_matrix={}".format(isf_matrix))

        # 修改_calc_schmidt_orthogonalization_tricky源码，得到指定算法的soe，并测试
        λ_ν_st = self.data_st.get_eigenvalue(self.ν_st)
        λ_ν = self.data_sn.get_eigenvalue(self.ν)
        e_value = λ_ν - λ_ν_st
        e_vectors = self.isf_func._calc_eigenspace_by_eigenvalue(isf_matrix, e_value)
        τ_max = len(e_vectors)

        logger.debug("e_vectors={}, τ_max={}".format(e_vectors, τ_max))

        soe_vectors = self.isf_func._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if τ_max > 1 \
            else [self.isf_func._calc_orthogonalization_vector(e_vectors[0])]

        soe_square = [Matrix([sp.sign(j) * j**2 for j in i]) for i in soe_vectors]
        logger.debug("soe_square={}".format(soe_square))
        '''
        小结：
        1，[3, 2, 1] * [3, 2, 1] = [3, 2, 1]的第一次尝试宣告失败
        细节：
        x1, x2, x3, x4, x5 = sp.symbols("x1 x2 x3 x4 x5")
        s增广矩阵：
        a = Matrix([[-sqrt(2)/4 - sqrt(2)/8, sqrt(30)/60 - -7*sqrt(30)/120, sqrt(6)/6 - sqrt(6)/6, 0 - 0,                   -sqrt(5)/5 - -sqrt(5)/5,    0], 
                    [sqrt(3)/4 - 0,          sqrt(5)/20 - sqrt(5)/5,        -Ra(1)/4 - -Ra(1)/4,   -sqrt(2)/4 - -sqrt(2)/4, -sqrt(30)/10 - sqrt(30)/10, 0],
                    [-sqrt(2)/8 - sqrt(2)/2, sqrt(30)/8 - 0,                0 - 0,                 0 - 0,                   0 - 0,                      0]])
        解：FiniteSet((sqrt(15)*x3/5, x2, x3, x4, 0))  # sp.linsolve(a, [x1, x2, x3, x4, x5])
        回拼的结果：
        Matrix([
            [        -3*sqrt(10)/20 - 2/5 - sqrt(5)/20],
            [                 -sqrt(30)/30 + sqrt(6)/6],
            [             -sqrt(2)/4 - 1/4 + sqrt(5)/5],
            [                 -sqrt(30)/30 + sqrt(6)/6],
            [-sqrt(15)/30 + 2*sqrt(3)/15 + sqrt(30)/10],
            [                 2*sqrt(2)/5 + sqrt(10)/5],
            [                              sqrt(30)/10],
            [             -sqrt(2)/4 - 1/4 + sqrt(5)/5],
            [                              sqrt(30)/10],
            [            -3*sqrt(5)/20 + 3*sqrt(10)/20]])
        符合S，但与书对不上（而且太复杂）
        a增广矩阵：
        b = Matrix([[-sqrt(2)/4 + sqrt(2)/8, sqrt(30)/60 + -7*sqrt(30)/120, sqrt(6)/6 + sqrt(6)/6, 0 + 0,                   -sqrt(5)/5 + -sqrt(5)/5,    0], 
                    [sqrt(3)/4 + 0,          sqrt(5)/20 + sqrt(5)/5,        -Ra(1)/4 + -Ra(1)/4,   -sqrt(2)/4 + -sqrt(2)/4, -sqrt(30)/10 + sqrt(30)/10, 0],
                    [-sqrt(2)/8 + sqrt(2)/2, sqrt(30)/8 + 0,                0 + 0,                 0 + 0,                   0 + 0,                      0]])
        解：FiniteSet((-sqrt(15)*x2/3, x2, sqrt(30)*x5/5, -sqrt(15)*x5/5, x5))  # sp.linsolve(b, [x1, x2, x3, x4, x5])
        
        2，[3, 2, 1] * [3, 2, 1] = [4, 2]的结果就是书中结果
        
        3，[3, 2, 1] * [3, 2, 1] = [4, 1, 1]的结果中，有1/4是直接命中的，其他需要重新组合
        细节见test_isf_321_321_411_41_τ4_vector
        '''

    # @pytest.mark.skip("pass")
    def test_all_sym_isf(self):
        """注意，没有assert，要看log"""
        # all_sym_key_list = ['σμν', 'μσν', 'νμσ', 'σνμ', 'νσμ', 'μνσ',
        #                     'σ~μ~ν', 'μ~σ~ν', 'ν~μ~σ', 'σ~ν~μ', 'ν~σ~μ', 'μ~ν~σ',
        #                     'σ~μν~', 'μ~σν~', 'ν~μσ~', 'σ~νμ~', 'ν~σμ~', 'μ~νσ~',
        #                     'σμ~ν~', 'μσ~ν~', 'νμ~σ~', 'σν~μ~', 'νσ~μ~', 'μν~σ~']
        all_sym_key_list = ['σμν', 'μσν', 'νμσ', 'σνμ', 'νσμ', 'μνσ']
        current_sym_key = self.sym_key
        for sym_key in all_sym_key_list:
            self.sym_key = sym_key
            logger.info("@@@@" * 8)
            logger.info("@@@@" * 8)
            logger.info("self.sym_key = {}".format(self.sym_key))
            logger.info("@@@@" * 4)
            logger.info("@@@@" * 4)
            self.test_make_sym_isf()
        self.sym_key = current_sym_key

    @pytest.mark.skip("pass")
    def test_make_sym_isf(self):
        """
        用于辅助生成书中没有的isf（使用书中结果和code中的对称函数）
        注意，没有assert，要看log
        """
        meta_σμν = (self.σ, self.μ, self.ν)
        sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(self.sym_key)
        sym_σμν = tuple(meta_σμν[d] if k is False else
                        self.data_sn.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
        sym_σ, sym_μ, sym_ν = sym_σμν
        logger.debug("meta:{}*{}={} by {}\n"
                     "sym:{}*{}={}\n"
                     "sym_key={}"
                     "".format(self.σ, self.μ, self.ν, self.ν_st, sym_σ, sym_μ, sym_ν, self.sym_key))
        meta_rows = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        sym_list = []  # (sym_row, sym_ν_st)
        for single_meta_row in meta_rows:
            meta_row_σ, meta_row_μ, meta_row_τ = single_meta_row if len(single_meta_row) == 3 \
                else (*single_meta_row, None)
            meta_σμν_st = (meta_row_σ, meta_row_μ, self.ν_st)

            # logger.warning("meta_σμν_st={}".format(meta_σμν_st))

            sym_σμν_st = tuple(meta_σμν_st[d] if k is False else
                               self.data_st.get_tilde(meta_σμν_st[d]) for d, k in zip(sym_d3, sym_k4))

            # logger.warning("sym_σμν_st={}".format(sym_σμν_st))

            single_sym_row = (sym_σμν_st[0], sym_σμν_st[1]) if meta_row_τ is None \
                else (sym_σμν_st[0], sym_σμν_st[1], meta_row_τ)
            sym_list_element = (single_sym_row, sym_σμν_st[2])
            sym_list.append(sym_list_element)
        logger.debug("meta_rows={}, sym_list={}".format(meta_rows, sym_list))
        sym_σμν = (sym_σ, sym_μ, sym_ν)
        sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(self.sym_key)
        helper = Helper(self.st, self.sn)
        for i, this_τ_meta_vector in enumerate(self.vector_list):
            if len(self.vector_list) == 1:
                sym_τ = None
            else:
                sym_τ = i + 1
            this_τ_meta_isf_square = [sp.sign(i) * (i ** 2) for i in this_τ_meta_vector]
            logger.info("\n{}*{}={}τ{} by {} =={}==> {}*{}={}τ{}, rows type below:"
                        "this_τ_meta_isf_square={}"
                        "".format(self.σ, self.μ, self.ν, sym_τ, self.ν_st, self.sym_key, sym_σ, sym_μ, sym_ν, sym_τ,
                                  this_τ_meta_isf_square))
            logger_msg = "{:<5}, {:<10} {:<25} {:<10}: " \
                         "{:<10} {:<25} {:<10} ?= {:<15}\n" \
                         "".format("Index", "meta_ν'", "meta_row", "meta_isf",
                                   "sym_ν'", "sym_row", "sym_isf", "sym_isf_code")

            # logger_msg = "{:<5}, {:<10} {:<30} = {:<15} from {:<30} = {:<16}\n" \
            #              "".format("Index", "sym_ν_st", "sym_row", "sym_isf_element",
            #                        "meta_row", "meta_isf_element")
            for (meta_i, meta_row), s_element in zip(enumerate(meta_rows), sym_list):
                sym_row, sym_ν_st = s_element
                sym_isf_element = helper.calc_sym_isf_element(sym_row, sym_σμν, sym_τ, sym_ν_st, self.sym_key,
                                                              sym_d3, sym_k4,
                                                              self.data_sn, self.data_st,
                                                              this_τ_meta_isf_square[meta_i])
                logger_msg += "{:<5}, {:<10} {:<25} {:<10}: " \
                              "{:<10} {:<25} {:<10} ?= {:<15}\n" \
                              "".format(meta_i, str(self.ν_st), str(meta_row), str(this_τ_meta_isf_square[meta_i]),
                                        str(sym_ν_st), str(sym_row), str(sym_isf_element), "SEE BOOK")
                # logger_msg += "{:<5}, {:<10} {:<30} = {:<15} from {:<30} = {:<16}\n" \
                #               "".format(meta_i, str(sym_ν_st), str(sym_row), str(sym_isf_element),
                #                         str(meta_row), str(this_τ_meta_isf_square[meta_i]))
            logger.info("\n{}\n".format(logger_msg))
            '''
            小结：分别利用4-195b式和对称式子，使用书中的[321]*[321]=[321]首相[32]，得到的结果不吻合！
            或者我的code错了，或者书中的做法，依然欠考虑（只考虑A~A~B，未考虑A~A~A~）
            
            对称式子部分，自己和自己倒是很自洽
            '''

    @pytest.mark.skip("pass")
    def test_isf_321_321_42_41_τ3_vector(self):
        """
        这个的结果是通过书中42_321_321_all对称来的
        ν' = [4, 1]                     τ1(s)       τ2(s)       τ3(s)
        ([3, 2], [3, 2])                0           5/16        25/512
        ([3, 2], [3, 1, 1])             3/20        -3/32       -3/1280
        ([3, 2], [2, 2, 1])             0           0           225/512
        ([3, 1, 1], [3, 2])             3/20        -3/32       -3/1280
        ([3, 1, 1], [3, 1, 1])          2/5         0           9/640
        ([3, 1, 1], [2, 2, 1])          3/20        3/32        -3/1280
        ([2, 2, 1], [3, 2])             0           0           225/512
        ([2, 2, 1], [3, 1, 1])          3/20        3/32        -3/1280
        ([2, 2, 1], [2, 2, 1])          0           -5/16       25/512
        """
        # 准备对称的部分，测试时就可以注释掉了
        # self.make_sym_isf()
        pass
        self.test_isf_calc_main_vector()
        '''
        结论：通过先[3, 2, 1] * [3, 2, 1] = [4, 2](τ=3)能够得到与书一致的首分支，继而能得到所有一致的对称
        '''

    @pytest.mark.skip("pass")
    def test_isf_321_321_411_41_τ4_vector(self):
        """
        这个的结果是通过书中411_321_321_all对称来的
        ν' = [4, 1]                     τ1(s)       τ2(a)       τ3(a)       τ4(a)
        ([3, 2], [3, 2])                75/256      0           0           0
        ([3, 2], [3, 1, 1])             5/128       3/512       1/4         125/512
        ([3, 2], [2, 2, 1])             -3/256      125/256     0           -3/256
        ([3, 1, 1], [3, 2])             5/128       -3/512      -1/4        -125/512
        ([3, 1, 1], [3, 1, 1])          -15/64      0           0           0
        ([3, 1, 1], [2, 2, 1])          5/128       3/512       -1/4        125/512
        ([2, 2, 1], [3, 2])             -3/256      -125/256    0           3/256
        ([2, 2, 1], [3, 1, 1])          5/128       -3/512      1/4         -125/512
        ([2, 2, 1], [2, 2, 1])          75/256      0           0           0
        """
        # 准备对称的部分，测试时就可以注释掉了
        # self.make_sym_isf()
        pass
        self.test_isf_calc_main_vector()
        '''
        最小soe：
        τ1 = Matrix([[5*sqrt(3)/16], [ sqrt(10)/16], [ -sqrt(3)/16], [ sqrt(10)/16], [ -sqrt(15)/8],
                     [ sqrt(10)/16], [ -sqrt(3)/16], [ sqrt(10)/16], [5*sqrt(3)/16]])
        τ2 = Matrix([[         0], [         0], [         0], [         0], [         0],
                     [-sqrt(2)/2], [         0], [ sqrt(2)/2], [         0]]) * -1             
        τ3 = Matrix([[         0], [         0], [-sqrt(2)/2], [         0], [         0],
                     [         0], [ sqrt(2)/2], [         0], [         0]]) * -1              
        τ4 = Matrix([[         0], [-sqrt(2)/2], [         0], [ sqrt(2)/2], [         0],
                     [         0], [         0], [         0], [         0]]) * -1  
        τ0 = [0] * len(τ1)
        τ_list = [τ1, τ2, τ3, τ4, τ0]
        x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
        s增广矩阵：(纯4.14)
        s = Matrix([[(i[0] - i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in τ_list]])  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
        解：FiniteSet((x1, -x4, 0, x4))  # sp.linsolve(s, [x1, x2, x3, x4])
        a增广矩阵：
        a = Matrix([[(i[0] + i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in τ_list]])  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
        解：FiniteSet((0, x4, x3, x4))  # sp.linsolve(a, [x1, x2, x3, x4])
        回拼的结果：
        s_1 = τ1
        a_1 = τ2 + τ4
        第一次重新正交归一化：
        sa_1_list = [s_1, a_1, τ3, τ4]
        sa_1_soe = sp.GramSchmidt(sa_1_list, True)
        [Matrix([[5*sqrt(3)/16],[ sqrt(10)/16],[ -sqrt(3)/16],[ sqrt(10)/16],[ -sqrt(15)/8],
                 [ sqrt(10)/16],[ -sqrt(3)/16],[ sqrt(10)/16],[5*sqrt(3)/16]]),
         Matrix([[   0],[ 1/2],[   0],[-1/2],[   0],[ 1/2],[   0],[-1/2],[   0]]),
         Matrix([[         0],[         0],[ sqrt(2)/2],[         0],[         0],
                 [         0],[-sqrt(2)/2],[         0],[         0]]),
         Matrix([[   0],[ 1/2],[   0],[-1/2],[   0],[-1/2],[   0],[ 1/2],[   0]])]
        平方：
        [75/256, 5/128, -3/256, 5/128, -15/64, 5/128, -3/256, 5/128, 75/256]
        [0, 1/4, 0, -1/4, 0, 1/4, 0, -1/4, 0]
        [0, 0, 1/2, 0, 0, 0, -1/2, 0, 0]
        [0, 1/4, 0, -1/4, 0, -1/4, 0, 1/4, 0]
        
        顺便看一下新一轮τ_list：
        new_τ_list = sa_1_soe[2::] + [τ0]
        s = Matrix([[(i[0] - i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in new_τ_list]])  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
        FiniteSet((0, x4))  # sp.linsolve(s, [x3, x4])
        a = Matrix([[(i[0] + i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in new_τ_list]])  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
        FiniteSet((x3, 0))  # sp.linsolve(a, [x3, x4])
        '''

    @pytest.mark.skip("pass")
    def test_isf_321_321_321_32_τ5_vector(self):
        """
        这个的结果就在书中，附上书中的对称性：τ1(s), τ2(a), τ3(a), τ4(s), τ5(s)
        """
        pass
        self.test_isf_calc_main_vector()
        '''
        最小soe:
        τ1 = Matrix([[-sqrt(15) / 20],     [-sqrt(2) / 4],  [sqrt(3) / 4],    [sqrt(2) / 8],         [sqrt(5) / 20],
                     [sqrt(30) / 20],      [-sqrt(2) / 8],  [0],              [sqrt(2) / 2],         [0]]) * -1
        τ2 = Matrix([[-Ra(1) / 4],         [sqrt(30) / 60], [sqrt(5) / 20],   [-7 * sqrt(30) / 120], [sqrt(3) / 12],
                     [sqrt(2) / 4],        [sqrt(30) / 8],  [sqrt(5) / 5],    [0],                   [0]])  * -1
        τ3 = Matrix([[-sqrt(5) / 20],      [sqrt(6) / 6],   [-Ra(1) / 4],     [sqrt(6) / 6],         [-sqrt(15) / 30],
                     [sqrt(10) / 5],       [0],             [-Ra(1) / 4],     [0],              [-3 * sqrt(5)/20]]) * -1
        τ4 = Matrix([[-3 * sqrt(10) / 20], [0],             [-sqrt(2) / 4],   [0],                   [sqrt(30) / 10],
                     [0],                  [0],             [-sqrt(2) / 4],   [0],              [3 * sqrt(10)/20]]) * -1
        τ5 = Matrix([[0],                  [-sqrt(5) / 5],  [-sqrt(30) / 10], [sqrt(5) / 5],         [0],
                     [0],                  [0],             [sqrt(30) / 10],  [0],                   [0]]) * -1
        τ0 = [0] * len(τ1)
        τ_list = [τ1, τ2, τ3, τ4, τ5, τ0]
        x1, x2, x3, x4, x5 = sp.symbols("x1 x2 x3 x4 x5")
        s增广矩阵：(纯4.14)
        s = Matrix([[(i[0] - i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] - i[4]) for i in τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] - i[5]) for i in τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        解：FiniteSet((sqrt(3)*x3/4, sqrt(5)*x3/4, x3, 0, 0))  # sp.linsolve(s, [x1, x2, x3, x4, x5])
        a增广矩阵：
        a = Matrix([[(i[0] + i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] + i[4]) for i in τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] + i[5]) for i in τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        解：FiniteSet((sqrt(10)*x5/2, -sqrt(6)*x5/2, 0, 0, x5))  # sp.linsolve(a, [x1, x2, x3, x4, x5])
        回拼的结果：
        s_1 = sqrt(3)/4 * τ1 + sqrt(5)/4 * τ2 + τ3
        a_1 = sqrt(10)/2 * τ1 + -sqrt(6)/2 * τ2 + τ5
        第一次重新正交归一化：
        sa_1_list = [s_1, a_1, τ3, τ4, τ5]
        sa_1_soe = sp.GramSchmidt(sa_1_list, True)
        [Matrix([[sqrt(30)/20],[       -1/4],[          0],[       -1/4],[          0],
                 [-sqrt(15)/5],[       -1/4],[          0],[       -1/4],[sqrt(30)/20]]),
         Matrix([[   0],[ 1/2],[   0],[-1/2],[   0],[   0],[ 1/2],[   0],[-1/2],[   0]]),
         Matrix([[-sqrt(15)/20],[  -sqrt(2)/4],[   sqrt(3)/4],[  -sqrt(2)/4],[  sqrt(5)/10],
                 [           0],[   sqrt(2)/4],[   sqrt(3)/4],[   sqrt(2)/4],[ sqrt(15)/20]]),
         Matrix([[ 3*sqrt(10)/20],[             0],[     sqrt(2)/4],[             0],[  -sqrt(30)/10],
                 [             0],[             0],[     sqrt(2)/4],[             0],[-3*sqrt(10)/20]]),
         Matrix([[         0],[       1/4],[ sqrt(6)/4],[      -1/4],[         0],
                 [         0],[      -1/4],[-sqrt(6)/4],[       1/4],[         0]])]
        平方的就是书中结果啦：
        [3/40, -1/16, 0, -1/16, 0, -3/5, -1/16, 0, -1/16, 3/40]   # s
        [0, 1/4, 0, -1/4, 0, 0, 1/4, 0, -1/4, 0]                  # a
        [-3/80, -1/8, 3/16, -1/8, 1/20, 0, 1/8, 3/16, 1/8, 3/80]  # a
        [9/40, 0, 1/8, 0, -3/10, 0, 0, 1/8, 0, -9/40]             # s
        [0, 1/16, 3/8, -1/16, 0, 0, -1/16, -3/8, 1/16, 0]         # s
        
        顺便看一下新一轮τ_list：
        new_τ_list = sa_1_soe[2::] + [τ0]
        s = Matrix([[(i[0] - i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in new_τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] - i[4]) for i in new_τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] - i[5]) for i in new_τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        FiniteSet((0, 0, 0))  # sp.linsolve(s, [x1, x2, x3])
        a = Matrix([[(i[0] + i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in new_τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] + i[4]) for i in new_τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] + i[5]) for i in new_τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        FiniteSet((0, 0, 0))  # sp.linsolve(a, [x1, x2, x3])
        '''

    @pytest.mark.skip("pass")
    def test_new_method_for_isf_321_321_321_32_τ5_vector(self):
        """
        论文中只考虑了AAB(A~=A)的情况，但面对AAA(A~=A)，实则应该有三个约束条件：
        1，μσν，这个约束同书中一致
        2，σ~μ~ν，这个约束书中也有
        3，νμσ，它涉及到sym结果要和首分支推非首分支一致的问题，是以前未曾考虑的
        注意：24个ϵ的独立约束只需要三个
        """
        # # 首先拿到首分支的soe
        # meta_rows = self.isf_func.calc_row_indexes(self.data_σ_μ, self.ν_st, self.data_st)
        # flag, isf_matrix = \
        #     self.isf_func._calc_isf_matrix_4_193a(meta_rows, self.ν_st, self.data_sn, self.data_σ_μ, self.data_st)
        # logger.debug("meta_rows={}, isf_matrix={}".format(meta_rows, isf_matrix))
        #
        # # 修改_calc_schmidt_orthogonalization_tricky源码，得到指定算法的soe，并测试
        # λ_ν_st = self.data_st.get_eigenvalue(self.ν_st)
        # λ_ν = self.data_sn.get_eigenvalue(self.ν)
        # e_value = λ_ν - λ_ν_st
        # e_vectors = self.isf_func._calc_eigenspace_by_eigenvalue(isf_matrix, e_value)
        # τ_max = len(e_vectors)
        #
        # logger.debug("e_vectors={}, τ_max={}".format(e_vectors, τ_max))
        #
        # soe_vectors = self.isf_func._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if τ_max > 1 \
        #     else [self.isf_func._calc_orthogonalization_vector(e_vectors[0])]
        #
        # soe_square = [Matrix([sp.sign(j) * j ** 2 for j in i]) for i in soe_vectors]
        # logger.debug("soe_square={}".format(soe_square))

        # 准备

        # 从sa_1_soe开始
        τ1 = Matrix([[-sqrt(15) / 20], [-sqrt(2) / 4], [sqrt(3) / 4], [sqrt(2) / 8], [sqrt(5) / 20],
                     [sqrt(30) / 20], [-sqrt(2) / 8], [0], [sqrt(2) / 2], [0]]) * -1
        τ2 = Matrix([[-Ra(1) / 4], [sqrt(30) / 60], [sqrt(5) / 20], [-7 * sqrt(30) / 120], [sqrt(3) / 12],
                     [sqrt(2) / 4], [sqrt(30) / 8], [sqrt(5) / 5], [0], [0]]) * -1
        τ3 = Matrix([[-sqrt(5) / 20], [sqrt(6) / 6], [-Ra(1) / 4], [sqrt(6) / 6], [-sqrt(15) / 30],
                     [sqrt(10) / 5], [0], [-Ra(1) / 4], [0], [-3 * sqrt(5) / 20]]) * -1
        τ4 = Matrix([[-3 * sqrt(10) / 20], [0], [-sqrt(2) / 4], [0], [sqrt(30) / 10],
                     [0], [0], [-sqrt(2) / 4], [0], [3 * sqrt(10) / 20]]) * -1
        τ5 = Matrix([[0], [-sqrt(5) / 5], [-sqrt(30) / 10], [sqrt(5) / 5], [0],
                     [0], [0], [sqrt(30) / 10], [0], [0]]) * -1
        s_1 = sqrt(3) / 4 * τ1 + sqrt(5) / 4 * τ2 + τ3
        a_1 = sqrt(10) / 2 * τ1 + -sqrt(6) / 2 * τ2 + τ5
        sa_1_list = [s_1, a_1, τ3, τ4, τ5]
        sa_1_soe = sp.GramSchmidt(sa_1_list, True)
        s_soe = [sa_1_soe[0], sa_1_soe[3], sa_1_soe[4]]
        a_soe = [sa_1_soe[1], sa_1_soe[2]]

        xa, xb, xc, xd, xe = sp.symbols("xa xb xc xd xe")
        # 下面只看s_line
        # 对于 ==νμσ==> 其对应是：
        # [3, 2] ([3, 1, 1], [3, 2])       ==νμσ==>: [3, 1, 1] ([3, 2], [3, 2])
        # [3, 2] ([3, 1, 1], [3, 1, 1], 1) ==νμσ==>: [3, 1, 1] ([3, 2], [3, 1, 1], 1)
        # [3, 2] ([3, 1, 1], [3, 1, 1], 2) ==νμσ==>: [3, 1, 1] ([3, 2], [3, 1, 1], 2)
        # [3, 2] ([3, 1, 1], [2, 2, 1])    ==νμσ==>: [3, 1, 1] ([3, 2], [2, 2, 1])
        # τ1   , τ4   , τ5               τ1   , τ4  , τ5
        # -1/16, 0    , 1/8    ==νμσ==>: -5/96, 0   , -5/48
        # 0    , -3/10, -1/20  ==νμσ==>: 0    , -1/4, 1/24
        # -3/5 , 0    , 0      ==νμσ==>: -1/2 , 0   , 0
        # -1/16, 0    , -1/8   ==νμσ==>: -5/96, 0   , 5/48
        # 因为它对称后还可以线性组合，所以，先对称，后组合
        # 简单可知，xa = -2 * xc

        # 对于4-195c，上面四个的结果是：
        # [3, 1, 1] ([3, 2], [3, 2])        : r1 = 0
        # [3, 1, 1] ([3, 2], [3, 1, 1], 1)  : r2
        # [3, 1, 1] ([3, 2], [3, 1, 1], 2)  : r3
        # [3, 1, 1] ([3, 2], [2, 2, 1])     : r4
        # r2 = 4/3*(-17*xa/160 - sqrt(6)*(sqrt(2)*xb/4 - sqrt(3)*xc/4)/8
        #           + 3*sqrt(10)*(-sqrt(30)*xb/10 - sqrt(5)*xc/10)/40
        #           + sqrt(30)*(sqrt(30)*xa/20 - 3*sqrt(10)*xb/20 - sqrt(15)*xc/20)/32
        #           - 7*sqrt(30)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/160)**2
        #      *sign(-17*xa/160 - sqrt(6)*(sqrt(2)*xb/4 - sqrt(3)*xc/4)/8
        #            + 3*sqrt(10)*(-sqrt(30)*xb/10 - sqrt(5)*xc/10)/40
        #            + sqrt(30)*(sqrt(30)*xa/20 - 3*sqrt(10)*xb/20 - sqrt(15)*xc/20)/32
        #            - 7*sqrt(30)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/160)
        # r3 = 4/3*(-3*sqrt(6)*xa/40 + 3*sqrt(2)*xb/32 - 3*sqrt(3)*xc/32
        #           + 3*sqrt(6)*(-xa/4 - sqrt(2)*xc/4)/16 - sqrt(6)*(-xa/4 + sqrt(2)*xc/4)/16
        #           + sqrt(15)*(-sqrt(30)*xb/10 - sqrt(5)*xc/10)/40
        #           - 3*sqrt(5)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/40)**2
        #      *sign(-3*sqrt(6)*xa/40 + 3*sqrt(2)*xb/32 - 3*sqrt(3)*xc/32
        #            + 3*sqrt(6)*(-xa/4 - sqrt(2)*xc/4)/16 - sqrt(6)*(-xa/4 + sqrt(2)*xc/4)/16
        #            + sqrt(15)*(-sqrt(30)*xb/10 - sqrt(5)*xc/10)/40
        #            - 3*sqrt(5)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/40)
        # r4 = 4/3*(-3*sqrt(10)*xa/40 + sqrt(10)*(-xa/4 - sqrt(2)*xc/4)/16
        #           + sqrt(10)*(-xa/4 + sqrt(2)*xc/4)/16
        #           - sqrt(3)*(sqrt(30)*xa/20 - 3*sqrt(10)*xb/20 - sqrt(15)*xc/20)/16
        #           - sqrt(3)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/16)**2
        #      *sign(-3*sqrt(10)*xa/40 + sqrt(10)*(-xa/4 - sqrt(2)*xc/4)/16
        #            + sqrt(10)*(-xa/4 + sqrt(2)*xc/4)/16
        #            - sqrt(3)*(sqrt(30)*xa/20 - 3*sqrt(10)*xb/20 - sqrt(15)*xc/20)/16
        #            - sqrt(3)*(sqrt(30)*xa/20 + 3*sqrt(10)*xb/20 + sqrt(15)*xc/20)/16)

        pass
        '''
        最小soe:
        τ1 = Matrix([[-sqrt(15) / 20],     [-sqrt(2) / 4],  [sqrt(3) / 4],    [sqrt(2) / 8],         [sqrt(5) / 20],
                     [sqrt(30) / 20],      [-sqrt(2) / 8],  [0],              [sqrt(2) / 2],         [0]]) * -1
        τ2 = Matrix([[-Ra(1) / 4],         [sqrt(30) / 60], [sqrt(5) / 20],   [-7 * sqrt(30) / 120], [sqrt(3) / 12],
                     [sqrt(2) / 4],        [sqrt(30) / 8],  [sqrt(5) / 5],    [0],                   [0]])  * -1
        τ3 = Matrix([[-sqrt(5) / 20],      [sqrt(6) / 6],   [-Ra(1) / 4],     [sqrt(6) / 6],         [-sqrt(15) / 30],
                     [sqrt(10) / 5],       [0],             [-Ra(1) / 4],     [0],              [-3 * sqrt(5)/20]]) * -1
        τ4 = Matrix([[-3 * sqrt(10) / 20], [0],             [-sqrt(2) / 4],   [0],                   [sqrt(30) / 10],
                     [0],                  [0],             [-sqrt(2) / 4],   [0],              [3 * sqrt(10)/20]]) * -1
        τ5 = Matrix([[0],                  [-sqrt(5) / 5],  [-sqrt(30) / 10], [sqrt(5) / 5],         [0],
                     [0],                  [0],             [sqrt(30) / 10],  [0],                   [0]]) * -1
        τ0 = [0] * len(τ1)
        τ_list = [τ1, τ2, τ3, τ4, τ5, τ0]
        x1, x2, x3, x4, x5 = sp.symbols("x1 x2 x3 x4 x5")
        s增广矩阵：(纯4.14)
        s = Matrix([[(i[0] - i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] - i[4]) for i in τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] - i[5]) for i in τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        解：FiniteSet((sqrt(3)*x3/4, sqrt(5)*x3/4, x3, 0, 0))  # sp.linsolve(s, [x1, x2, x3, x4, x5])
        a增广矩阵：
        a = Matrix([[(i[0] + i[-1]) for i in τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] + i[4]) for i in τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] + i[5]) for i in τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        解：FiniteSet((sqrt(10)*x5/2, -sqrt(6)*x5/2, 0, 0, x5))  # sp.linsolve(a, [x1, x2, x3, x4, x5])
        回拼的结果：
        s_1 = sqrt(3)/4 * τ1 + sqrt(5)/4 * τ2 + τ3
        a_1 = sqrt(10)/2 * τ1 + -sqrt(6)/2 * τ2 + τ5
        第一次重新正交归一化：
        sa_1_list = [s_1, a_1, τ3, τ4, τ5]
        sa_1_soe = sp.GramSchmidt(sa_1_list, True)
        [Matrix([[sqrt(30)/20],[       -1/4],[          0],[       -1/4],[          0],
                 [-sqrt(15)/5],[       -1/4],[          0],[       -1/4],[sqrt(30)/20]]),
         Matrix([[   0],[ 1/2],[   0],[-1/2],[   0],[   0],[ 1/2],[   0],[-1/2],[   0]]),
         Matrix([[-sqrt(15)/20],[  -sqrt(2)/4],[   sqrt(3)/4],[  -sqrt(2)/4],[  sqrt(5)/10],
                 [           0],[   sqrt(2)/4],[   sqrt(3)/4],[   sqrt(2)/4],[ sqrt(15)/20]]),
         Matrix([[ 3*sqrt(10)/20],[             0],[     sqrt(2)/4],[             0],[  -sqrt(30)/10],
                 [             0],[             0],[     sqrt(2)/4],[             0],[-3*sqrt(10)/20]]),
         Matrix([[         0],[       1/4],[ sqrt(6)/4],[      -1/4],[         0],
                 [         0],[      -1/4],[-sqrt(6)/4],[       1/4],[         0]])]
        平方的就是书中结果啦：
        [3/40, -1/16, 0, -1/16, 0, -3/5, -1/16, 0, -1/16, 3/40]   # s
        [0, 1/4, 0, -1/4, 0, 0, 1/4, 0, -1/4, 0]                  # a
        [-3/80, -1/8, 3/16, -1/8, 1/20, 0, 1/8, 3/16, 1/8, 3/80]  # a
        [9/40, 0, 1/8, 0, -3/10, 0, 0, 1/8, 0, -9/40]             # s
        [0, 1/16, 3/8, -1/16, 0, 0, -1/16, -3/8, 1/16, 0]         # s
        
        顺便看一下新一轮τ_list：
        new_τ_list = sa_1_soe[2::] + [τ0]
        s = Matrix([[(i[0] - i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] - i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] - i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] - i[-4]) for i in new_τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] - i[4]) for i in new_τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] - i[5]) for i in new_τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        FiniteSet((0, 0, 0))  # sp.linsolve(s, [x1, x2, x3])
        a = Matrix([[(i[0] + i[-1]) for i in new_τ_list],  # ([3, 2], [3, 2]) <-> ([2, 2, 1], [2, 2, 1])
                    [(i[1] + i[-2]) for i in new_τ_list],  # ([3, 2], [3, 1, 1]) <-> ([2, 2, 1], [3, 1, 1])
                    [(i[2] + i[-3]) for i in new_τ_list],  # ([3, 2], [2, 2, 1]) <-> ([2, 2, 1], [3, 2])
                    [(i[3] + i[-4]) for i in new_τ_list],  # ([3, 1, 1], [3, 2]) <-> ([3, 1, 1], [2, 2, 1])
                    [(i[4] + i[4]) for i in new_τ_list],  # ([3, 1, 1], [3, 1, 1], 1) <-> ([3, 1, 1], [3, 1, 1], 1) 
                    [(i[5] + i[5]) for i in new_τ_list]])  # ([3, 1, 1], [3, 1, 1], 2) <-> ([3, 1, 1], [3, 1, 1], 2) 
        FiniteSet((0, 0, 0))  # sp.linsolve(a, [x1, x2, x3])
        
        再将上面全改造为s/a的分组，再做一次线性组合，使其满足式子4-195c
        s_soe = [sa_1_soe[0], sa_1_soe[3], sa_1_soe[4]]
        a_soe = [sa_1_soe[1], sa_1_soe[2]]
        '''

    @pytest.mark.skip("pass")
    def test_ϵ_31_211_31_νμσ(self):
        """顺便查一下ϵ，不作为保留test"""
        '''
        代码：
        code_dict = {'σμν': 1,'σ~μ~ν': -1,'μσν': 1,'μ~σ~ν': -1,'σ~μν~': 1,
                     'σμ~ν~': 1,'μ~σν~': -1,'μσ~ν~': -1,'σνμ': -1,'σ~ν~μ': -1,
                     'νσμ': 1,'ν~σ~μ': 1,'σ~νμ~': 1,'σν~μ~': -1,'ν~σμ~': -1,
                     'νσ~μ~': 1,'νμσ': -1,'ν~μ~σ': 1,'μνσ': -1,'μ~ν~σ': 1,
                     'ν~μσ~': -1,'νμ~σ~': -1,'μ~νσ~': 1,'μν~σ~': 1}
        code_flags = {'σμν': (2, 1, 1),'σ~μ~ν': (3, 2, 1),'μσν': (2, 1, 1),'μ~σ~ν': (3, 2, 1),'σ~μν~': (2, 3, 3),
                      'σμ~ν~': (1, 2, 3),'μ~σν~': (2, 3, 3),'μσ~ν~': (1, 2, 3),'σνμ': (1, 1, 2),'σ~ν~μ': (2, 1, 1),
                      'νσμ': (2, 1, 1),'ν~σ~μ': (1, 1, 2),'σ~νμ~': (3, 3, 2),'σν~μ~': (2, 3, 3),'ν~σμ~': (2, 3, 3),
                      'νσ~μ~': (3, 3, 2),'νμσ': (1, 1, 2),'ν~μ~σ': (1, 2, 3),'μνσ': (1, 1, 2),'μ~ν~σ': (1, 2, 3),
                      'ν~μσ~': (3, 3, 2),'νμ~σ~': (3, 2, 1),'μ~νσ~': (3, 3, 2),'μν~σ~': (3, 2, 1)}
        '''
        meta_ϵ_dict, meta_ϵ_flags = self.ϵ_func.calc_meta_ϵ_dict(self.data_σ_μ, self.ν, None, self.data_sn)
        meta_data_σ_μ, meta_ν = self.data_σ_μ, self.ν
        meta_yd_σμν = (meta_data_σ_μ.σ, meta_data_σ_μ.μ, meta_ν)
        meta_Λ_list_list = [self.data_sn.get_phase_factor_list(yd) for yd in meta_yd_σμν]
        meta_h_list = [self.data_sn.get_yt_num(yd) for yd in meta_yd_σμν]
        logger.warning("meta_yd_σμν={}".format(meta_yd_σμν))
        logger.warning("meta_ϵ_dict = {}".format(meta_ϵ_dict))
        logger.warning("meta_ϵ_flags = {}".format(meta_ϵ_flags))
        '''
        直接计算[3, 1] * [2, 1, 1] = [3, 1]:
        ϵ_dict = {'σμν': 1, 'σ~μ~ν': -1, 'μσν': 1, 'μ~σ~ν': -1, 'σ~μν~': 1, 'σμ~ν~': 1, 'μ~σν~': -1, 'μσ~ν~': -1, 
                  'σνμ': -1, 'σ~ν~μ': -1, 'νσμ': 1, 'ν~σ~μ': 1, 'σ~νμ~': 1, 'σν~μ~': -1, 'ν~σμ~': -1, 'νσ~μ~': 1, 
                  'νμσ': -1, 'ν~μ~σ': 1, 'μνσ': -1, 'μ~ν~σ': 1, 'ν~μσ~': -1, 'νμ~σ~': -1, 'μ~νσ~': 1, 'μν~σ~': 1}
        ϵ_flags = {'σμν': (2, 1, 1), 'σ~μ~ν': (3, 2, 1), 'μσν': (2, 1, 1), 'μ~σ~ν': (3, 2, 1), 'σ~μν~': (2, 3, 3), 
                   'σμ~ν~': (1, 2, 3), 'μ~σν~': (2, 3, 3), 'μσ~ν~': (1, 2, 3), 'σνμ': (1, 1, 2), 'σ~ν~μ': (2, 1, 1), 
                   'νσμ': (2, 1, 1), 'ν~σ~μ': (1, 1, 2), 'σ~νμ~': (3, 3, 2), 'σν~μ~': (2, 3, 3), 'ν~σμ~': (2, 3, 3), 
                   'νσ~μ~': (3, 3, 2), 'νμσ': (1, 1, 2), 'ν~μ~σ': (1, 2, 3), 'μνσ': (1, 1, 2), 'μ~ν~σ': (1, 2, 3), 
                   'ν~μσ~': (3, 3, 2), 'νμ~σ~': (3, 2, 1), 'μ~νσ~': (3, 3, 2), 'μν~σ~': (3, 2, 1)}
        '''
        # tmp_sym_ϵ_dict = None
        # tmp_sym_ϵ_flags = None
        # for sym_mode, sym_mode_d3, sym_mode_k4 in \
        #         chain(self.ϵ_func.σμν_0, self.ϵ_func.σμν_1, self.ϵ_func.σμν_2,
        #               self.ϵ_func.σμν_3, self.ϵ_func.σμν_4, self.ϵ_func.σμν_5):
        #     sym_yd_σμν = tuple(meta_yd_σμν[d] if k is False else
        #                        self.data_sn.get_tilde(meta_yd_σμν[d]) for d, k in zip(sym_mode_d3, sym_mode_k4))
        #     if sym_yd_σμν != ([3, 1], [2, 1, 1], [3, 1]):
        #         continue
        #     logger.debug("\n ===================")
        #     logger.warning("sym_mode={}, sym_mode_d3={}, sym_mode_k4={}".format(sym_mode, sym_mode_d3, sym_mode_k4))
        #     sym_Λ_list_list = [self.data_sn.get_phase_factor_list(yd) for yd in sym_yd_σμν]
        #     sym_ϵ_dict, sym_ϵ_flags = \
        #         self.ϵ_func._calc_single_sym_mode_ϵ_by_meta(sym_mode, sym_mode_d3, sym_mode_k4, meta_ϵ_dict, meta_ϵ_flags,
        #                                                     meta_Λ_list_list, sym_Λ_list_list, meta_h_list, self.data_sn)
        #     logger.warning("sym_yd_σμν={}".format(sym_yd_σμν))
        #     if tmp_sym_ϵ_dict is None and tmp_sym_ϵ_flags is None:
        #         tmp_sym_ϵ_dict = sym_ϵ_dict
        #         tmp_sym_ϵ_flags = sym_ϵ_flags
        #         logger.warning("tmp_sym_ϵ_dict = {}".format(tmp_sym_ϵ_dict))
        #         logger.warning("tmp_sym_ϵ_flags = {}".format(tmp_sym_ϵ_flags))
        #     else:
        #         if tmp_sym_ϵ_dict == sym_ϵ_dict and tmp_sym_ϵ_flags == sym_ϵ_flags:
        #             continue
        #         else:
        #             logger.warning("sym_ϵ_dict = {}".format(sym_ϵ_dict))
        #             logger.warning("sym_ϵ_flags = {}".format(sym_ϵ_flags))
        '''
        先计算meta[3, 1] * [3, 1] = [2, 1, 1] ，再对称成[3, 1] * [2, 1, 1] = [3, 1]:
        a = {'σμν': 1, 'σ~μ~ν': -1, 'μσν': 1, 'μ~σ~ν': -1, 'σ~μν~': 1, 'σμ~ν~': 1, 'μ~σν~': -1, 'μσ~ν~': -1, 
             'σνμ': -1, 'σ~ν~μ': -1, 'νσμ': 1, 'ν~σ~μ': 1, 'σ~νμ~': 1, 'σν~μ~': -1, 'ν~σμ~': -1, 'νσ~μ~': 1, 
             'νμσ': -1, 'ν~μ~σ': 1, 'μνσ': -1, 'μ~ν~σ': 1, 'ν~μσ~': -1, 'νμ~σ~': -1, 'μ~νσ~': 1, 'μν~σ~': 1}
        b = {'σμν': (2, 1, 1), 'σ~μ~ν': (3, 2, 1), 'μσν': (2, 1, 1), 'μ~σ~ν': (3, 2, 1), 'σ~μν~': (2, 3, 3), 
             'σμ~ν~': (1, 2, 3), 'μ~σν~': (2, 3, 3), 'μσ~ν~': (1, 2, 3), 'σνμ': (1, 1, 2), 'σ~ν~μ': (2, 1, 1), 
             'νσμ': (2, 1, 1), 'ν~σ~μ': (1, 1, 2), 'σ~νμ~': (3, 3, 2), 'σν~μ~': (2, 3, 3), 'ν~σμ~': (2, 3, 3), 
             'νσ~μ~': (3, 3, 2), 'νμσ': (1, 1, 2), 'ν~μ~σ': (1, 2, 3), 'μνσ': (1, 1, 2), 'μ~ν~σ': (1, 2, 3), 
             'ν~μσ~': (3, 3, 2), 'νμ~σ~': (3, 2, 1), 'μ~νσ~': (3, 3, 2), 'μν~σ~': (3, 2, 1)}
        它和直接计算的结果是一致的
        '''


class Helper(object):
    """这里放一些辅助函数"""

    def __init__(self, st, sn):
        self.st = st
        self.sn = sn

    def calc_sym_isf_element(self, s_row, sym_σμν, sym_τ, sym_ν_st, sym_key, sym_d3, sym_k4,
                             data_sn, data_st, single_meta_isf_square):
        """改自isf_and_cgc.py:_calc_sym_isf_element
        只改变必要的变量获取途径，不改变算法"""
        # 0, 数据准备
        sym_σ_st, sym_μ_st, sym_τ_st = s_row if len(s_row) == 3 else (s_row[0], s_row[1], None)
        meta_τ_st = sym_τ_st
        sym_σμν_st = (sym_σ_st, sym_μ_st, sym_ν_st)
        # d3和k4对于σμν，σ'μ'ν'应该是一致的，所以用它也能反推meta'
        meta_σμν = [None, None, None]
        meta_σμν_st = [None, None, None]
        for single_s_σμν, single_s_σμν_st, s_d, s_k in zip(sym_σμν, sym_σμν_st, sym_d3, sym_k4):
            single_meta_σμν = single_s_σμν if s_k is False else data_sn.get_tilde(single_s_σμν)
            meta_σμν[s_d] = single_meta_σμν
            single_meta_σμν_st = single_s_σμν_st if s_k is False else data_st.get_tilde(single_s_σμν_st)
            meta_σμν_st[s_d] = single_meta_σμν_st
        meta_τ = sym_τ

        # 1, h_ν'*h_μ / h_ν*h_μ'
        h_meta_ν = data_sn.get_yt_num(meta_σμν[-1])
        h_meta_ν_st = data_st.get_yt_num(meta_σμν_st[-1])
        h_sym_ν = data_sn.get_yt_num(sym_σμν[-1])
        h_sym_ν_st = data_st.get_yt_num(sym_σμν_st[-1])
        h_square = sp.Rational(h_meta_ν_st * h_sym_ν, h_meta_ν * h_sym_ν_st)

        # 2, ϵ
        # ϵ_dict, _ = meta_data_σ_μ.get_ϵ_dict_and_flags(meta_σμν[-1], meta_τ)
        # ϵ_dict也要load获取
        _, rst_dict = load_ϵ(*(self.sn, *meta_σμν, meta_τ), is_with_flags=True)
        ϵ_dict = rst_dict.get("data", {})
        _, rst_st_dict = load_ϵ(*(self.st, *meta_σμν_st, meta_τ_st), is_with_flags=True)
        ϵ_st_dict, ϵ_st_flags = rst_st_dict.get("data", {}), rst_st_dict.get("flags", {})
        ϵ_4 = ϵ_st_dict["σμν"] * ϵ_dict["σμν"] * ϵ_st_dict[sym_key] * ϵ_dict[sym_key]

        # 3, Λ
        meta_st_Λ_list_list = [data_st.get_phase_factor_list(yd) for yd in meta_σμν_st]
        meta_Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in meta_σμν]
        meta_m_st_list = ϵ_st_flags[sym_key]  # 其实取1也可以，m'是任意取的
        meta_m_list = [data_sn.quick_calc_m(m_st, yd_st, yd)
                       for m_st, yd_st, yd in zip(meta_m_st_list, meta_σμν_st, meta_σμν)]
        ΛΛ = 1
        for d, k in zip(sym_d3, sym_k4):
            if k is True:
                ΛΛ *= meta_Λ_list_list[d][meta_m_list[d] - 1] * meta_st_Λ_list_list[d][meta_m_st_list[d] - 1]

        # 4, ISF_σσ'μμ'νν'
        # 跳过读取，直接传入single_meta_isf_square

        # ISF_ν~ν'~σσ'μ~μ'~

        # logger.warning("\ns_row={}, sym_key={}".format(s_row, sym_key))
        # logger.warning("meta_σμν={}, ϵ_dict={}\n"
        #                "meta_σμν_st={}, ϵ_st_dict={}".format(meta_σμν, ϵ_dict, meta_σμν_st, ϵ_st_dict))
        # logger.warning("h_square={}, ϵ_4={}, ΛΛ={}".format(h_square, ϵ_4, ΛΛ))

        single_sym_isf_square = h_square * ϵ_4 * ΛΛ * single_meta_isf_square

        return single_sym_isf_square


'''
σ=[5, 1, 1], μ=[4, 2, 1] ν_st=[3, 1, 1, 1] isf_matrix的生成细节 它们的S5 isf都通过了测试
row_index_list = [([5, 1], [4, 1, 1]), ([5, 1], [3, 2, 1]), ([4, 1, 1], [4, 2]), ([4, 1, 1], [4, 1, 1]), 
                  ([4, 1, 1], [3, 2, 1], 1), ([4, 1, 1], [3, 2, 1], 2)]
col_index_list = [([4, 1, 1, 1], 1), ([4, 1, 1, 1], 2), ([3, 2, 1, 1], 1), ([3, 2, 1, 1], 2), ([3, 2, 1, 1], 3),
                  [3, 1, 1, 1, 1]]
                

1, method1

isf_matrix_with_4_193a=Matrix(
[[1/3, -2*sqrt(2)/3, -sqrt(35)/5, -sqrt(7)/3, 0, -4*sqrt(35)/15], 
 [-2*sqrt(2)/3, -5/6, -sqrt(70)/5, -sqrt(14)/3, sqrt(14)/4, -sqrt(70)/60], 
 [-sqrt(35)/5, -sqrt(70)/5, -7/5, -sqrt(5)/5, -sqrt(5)/5, -1/5], 
 [-sqrt(7)/3, -sqrt(14)/3, -sqrt(5)/5, 1/3, -1, 29*sqrt(5)/75], 
 [0, sqrt(14)/4, -sqrt(5)/5, -1, 7/4, 3*sqrt(5)/20], 
 [-4*sqrt(35)/15, -sqrt(70)/60, -1/5, 29*sqrt(5)/75, 3*sqrt(5)/20, 109/60]])
！ 它的本征值不对，无解（ 29*sqrt(5)/75 这个矩阵元和下面的不同！）
 
isf_matrix_with_last_m=Matrix(
[[1/3, -2*sqrt(2)/3, -sqrt(35)/5, -sqrt(7)/3, 0, -4*sqrt(35)/15], 
 [-2*sqrt(2)/3, -5/6, -sqrt(70)/5, -sqrt(14)/3, sqrt(14)/4, -sqrt(70)/60], 
 [-sqrt(35)/5, -sqrt(70)/5, -7/5, -sqrt(5)/5, -sqrt(5)/5, -1/5], 
 [-sqrt(7)/3, -sqrt(14)/3, -sqrt(5)/5, 1/3, -1, sqrt(5)/15], 
 [0, sqrt(14)/4, -sqrt(5)/5, -1, 7/4, 3*sqrt(5)/20], 
 [-4*sqrt(35)/15, -sqrt(70)/60, -1/5, sqrt(5)/15, 3*sqrt(5)/20, 109/60]])
$ 它看起来是有效的（不敢说正确的，但起码本征值对上了，本征矢量也是实矢量）

isf_matrix_with_first_m=Matrix(
[[1/3, -2*sqrt(2)/3, -sqrt(35)/5, -sqrt(7)/3, 0, -8*sqrt(35)/45], 
 [-2*sqrt(2)/3, -5/6, -sqrt(70)/5, -sqrt(14)/3, sqrt(14)/4, -sqrt(70)/90], 
 [-sqrt(35)/5, -sqrt(70)/5, -7/5, -sqrt(5)/5, -sqrt(5)/5, -2/15], 
 [-sqrt(7)/3, -sqrt(14)/3, -sqrt(5)/5, 1/3, -1, 2*sqrt(5)/45], 
 [0, sqrt(14)/4, -sqrt(5)/5, -1, 7/4, sqrt(5)/10], 
 [-8*sqrt(35)/45, -sqrt(70)/90, -2/15, 2*sqrt(5)/45, sqrt(5)/10, 341/180]])
! 它的本征矢量不对，出现了虚数

$$$ 从isf_matrix_with_4_193a和isf_matrix_with_last_m看出，
    对应到row_index_list，就是left_st=([4, 1, 1], [4, 1, 1])，right_st=([4, 1, 1], [3, 2, 1], 2)时出了问题
    奇怪的是，如果right_st=([4, 1, 1], [3, 2, 1], 2)有问题，([4, 1, 1], [3, 2, 1], 1)居然没有问题？

####_calc_isf_matrix_4_193a####
h_ν_st = 10
ν_st_st=[3, 1, 1], h_ν_st_st=6, 贡献6*6/10 = 18/5
sub_isf_matrix=Matrix(
[[1/18, -2*sqrt(2)/9, -sqrt(35)/18, -sqrt(7)/18, 0, 0], 
 [-2*sqrt(2)/9, -5/36, -sqrt(70)/18, -sqrt(14)/18, -sqrt(14)/72, -sqrt(70)/72], 
 [-sqrt(35)/18, -sqrt(70)/18, -17/90, -11*sqrt(5)/90, sqrt(5)/90, -11/90], 
 [-sqrt(7)/18, -sqrt(14)/18, -11*sqrt(5)/90, 1/6, 1/18, sqrt(5)/18], 
 [0, -sqrt(14)/72, sqrt(5)/90, 1/18, 17/72, -sqrt(5)/120], 
 [0, -sqrt(70)/72, -11/90, sqrt(5)/18, -sqrt(5)/120, 73/360]])
 
ν_st_st=[2, 1, 1, 1], h_ν_st_st=4, 贡献4*6/10 = 12/5
sub_isf_matrix=Matrix(
[[1/18, sqrt(2)/18, 0, -sqrt(7)/18, 0, -sqrt(35)/9], 
 [sqrt(2)/18, -5/36, 0, -sqrt(14)/18, sqrt(14)/8, sqrt(70)/72], 
 [0, 0, -3/10, sqrt(5)/10, -sqrt(5)/10, 1/10], 
 [-sqrt(7)/18, -sqrt(14)/18, sqrt(5)/10, -1/9, -1/2, 7*sqrt(5)/90], 
 [0, sqrt(14)/8, -sqrt(5)/10, -1/2, 3/8, 3*sqrt(5)/40], 
 [-sqrt(35)/9, sqrt(70)/72, 1/10, 7*sqrt(5)/90, 3*sqrt(5)/40, 163/360]])
 
####_calc_isf_matrix_with_last_m####
i_n=(1, 7), sub_isf_matrix=Matrix(
[[1/18, -7*sqrt(2)/36, -sqrt(35)/20, -sqrt(7)/18, 0, -sqrt(35)/90], 
 [-7*sqrt(2)/36, -5/36, -sqrt(70)/20, -sqrt(14)/18, 0, -sqrt(70)/90], 
 [-sqrt(35)/20, -sqrt(70)/20, -1/5, -sqrt(5)/10, 0, -1/10], 
 [-sqrt(7)/18, -sqrt(14)/18, -sqrt(5)/10, 5/36, 0, -sqrt(5)/45], 
 [0, 0, 0, 0, 1/4, 0], 
 [-sqrt(35)/90, -sqrt(70)/90, -1/10, -sqrt(5)/45, 0, 41/180]])

i_n=(2, 7),(3, 7),(4, 7), sub_isf_matrix = same as i_n=(1, 7)

i_n=(5, 7), sub_isf_matrix=Matrix(
[[1/18, sqrt(2)/18, 0, -sqrt(7)/18, 0, -sqrt(35)/9], 
 [sqrt(2)/18, -5/36, 0, -sqrt(14)/18, sqrt(14)/8, sqrt(70)/72], 
 [0, 0, -3/10, sqrt(5)/10, -sqrt(5)/10, 1/10], 
 [-sqrt(7)/18, -sqrt(14)/18, sqrt(5)/10, -1/9, -1/2, 7*sqrt(5)/90], 
 [0, sqrt(14)/8, -sqrt(5)/10, -1/2, 3/8, 3*sqrt(5)/40], 
 [-sqrt(35)/9, sqrt(70)/72, 1/10, 7*sqrt(5)/90, 3*sqrt(5)/40, 163/360]])

i_n=(6, 7), sub_isf_matrix = same as i_n=(5, 7)

####_calc_isf_matrix_with_first_m####
i_n=(1, 7), sub_isf_matrix=Matrix(
[[1/18, 0, -sqrt(35)/90, -sqrt(7)/18, 0, -8*sqrt(35)/135], 
 [0, -5/36, -sqrt(70)/90, -sqrt(14)/18, 7*sqrt(14)/72, 7*sqrt(70)/1080], 
 [-sqrt(35)/90, -sqrt(70)/90, -5/18, sqrt(5)/18, -7*sqrt(5)/90, 1/30], 
 [-sqrt(7)/18, -sqrt(14)/18, sqrt(5)/18, -1/18, -7/18, sqrt(5)/30], 
 [0, 7*sqrt(14)/72, -7*sqrt(5)/90, -7/18, 25/72, 49*sqrt(5)/1080], 
 [-8*sqrt(35)/135, 7*sqrt(70)/1080, 1/30, sqrt(5)/30, 49*sqrt(5)/1080, 439/1080]])

i_n=(2, 7), (3, 7), sub_isf_matrix = sub_isf_matrix = same as i_n=(1, 7)

i_n=(4, 7), sub_isf_matrix=Matrix(
[[1/18, -2*sqrt(2)/9, -sqrt(35)/18, -sqrt(7)/18, 0, 0], 
 [-2*sqrt(2)/9, -5/36, -sqrt(70)/18, -sqrt(14)/18, -sqrt(14)/72, -sqrt(70)/120], 
 [-sqrt(35)/18, -sqrt(70)/18, -17/90, -11*sqrt(5)/90, sqrt(5)/90, -1/18], 
 [-sqrt(7)/18, -sqrt(14)/18, -11*sqrt(5)/90, 1/6, 1/18, -sqrt(5)/18], 
 [0, -sqrt(14)/72, sqrt(5)/90, 1/18, 17/72, -sqrt(5)/72], 
 [0, -sqrt(70)/120, -1/18, -sqrt(5)/18, -sqrt(5)/72, 17/72]])

i_n=(5, 7), sub_isf_matrix = sub_isf_matrix = same as i_n=(4, 7)

i_n=(6, 7), sub_isf_matrix=Matrix(
[[1/18, -2*sqrt(2)/9, -sqrt(35)/18, -sqrt(7)/18, 0, 0], 
 [-2*sqrt(2)/9, -5/36, -sqrt(70)/18, -sqrt(14)/18, -sqrt(14)/72, -sqrt(70)/72], 
 [-sqrt(35)/18, -sqrt(70)/18, -17/90, -11*sqrt(5)/90, sqrt(5)/90, -11/90], 
 [-sqrt(7)/18, -sqrt(14)/18, -11*sqrt(5)/90, 1/6, 1/18, sqrt(5)/18], 
 [0, -sqrt(14)/72, sqrt(5)/90, 1/18, 17/72, -sqrt(5)/120], 
 [0, -sqrt(70)/72, -11/90, sqrt(5)/18, -sqrt(5)/120, 73/360]])
 
$$$$ 小结 $$$$
看起来问题出现在_calc_isf_matrix_with_first_m的i_n=(6, 7), sub_isf_matrix中
'''
