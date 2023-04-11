# -*- coding:utf8 -*-
"""
验收:
所有isf的绝对相位（英文书4.19.2 Overall phase）
所有isf满足式4-195b（从首项推非首项）
所有isf的所有对称自洽
"""


import pytest
import sympy as sp
from sympy import Rational as Ra
from sympy import sqrt, Matrix
from itertools import product, combinations, combinations_with_replacement, chain, permutations
from conf.cgc_config import default_s_n, group_d3, group_k4
from conf.cgc_config import min_s_n_of_ϵ, min_s_n_of_isf
from core.young_diagrams import load_young_diagrams
from core.young_tableaux import load_young_table_num
from core.branching_laws import load_branching_law
from core.cg_series import load_cg_series
from core.symmetry_combination import calc_ϵ_map_dicts
from core.isf_and_cgc import load_isf, DataHelper, ISFHelper, ΣMDataHelper, load_ϵ
from utils.log import get_logger


logger = get_logger(__name__)


class BaseCls(object):

    def setup_class(self):
        # 它不改变数据库，所以不需要保护
        # 准备
        self.sn_min = max(min_s_n_of_ϵ, min_s_n_of_isf, 3)
        self.sn_max = 5
        self.data_st = None
        self.data_sn = None
        self.isf_func = None
        self.data_σ_μ = None
        self.ϵ_key2groups_dict, self.groups2ϵ_key_dict = calc_ϵ_map_dicts()
        self.all_sym_key_list = ['σμν', 'μσν', 'νμσ', 'σνμ', 'νσμ', 'μνσ',
                                 'σ~μ~ν', 'μ~σ~ν', 'ν~μ~σ', 'σ~ν~μ', 'ν~σ~μ', 'μ~ν~σ',
                                 'σ~μν~', 'μ~σν~', 'ν~μσ~', 'σ~νμ~', 'ν~σμ~', 'μ~νσ~',
                                 'σμ~ν~', 'μσ~ν~', 'νμ~σ~', 'σν~μ~', 'νσ~μ~', 'μν~σ~']
        self.preload_isf_dict = None
        self.isf_2_self_sym_isf_format = "\n{:<10} {:<25} {:<10} {:<5}: {:<10} {:<25} {:<10} != {:<15}"
        self.isf_2_self_sym_isf_head = self.isf_2_self_sym_isf_format.format(
            "meta_ν_st", "meta_row", "meta_isf", "τ", "sym_ν_st", "sym_row", "sym_isf", "sym_isf_code")

        # 防bug
        '''
        这个bug来自sympy1.0版本之后的类注册器。如果不手动建立一次非整数的Rational，让它的注册表里存在Rational，直接读pickle会报错。

        C, including its class ClassRegistry, has been deprecated since SymPy
        1.0. It will be last supported in SymPy version 1.0. Use direct
        imports from the defining module instead. See
        https://github.com/sympy/sympy/issues/9371 for more info.
        '''
        _ = Matrix([[Ra(1) / 2]])

    def teardown_class(self):
        pass

    def _prepare_sn(self, s_n):
        s_t = s_n - 1
        self.data_st = DataHelper(s_t, self.ϵ_key2groups_dict, self.groups2ϵ_key_dict)
        data_st_st_yt_num_list = []
        _, yd_st_st_list = load_young_diagrams(s_n - 1 - 1)
        for yd_st_st in yd_st_st_list:
            _, total_num = load_young_table_num(s_n - 1 - 1, yd_st_st)
            data_st_st_yt_num_list.append(total_num)
        self.data_st.add_more_data(data_st_st_yt_num_list, yd_st_st_list)

        self.data_sn = DataHelper(s_n, self.ϵ_key2groups_dict, self.groups2ϵ_key_dict)
        self.data_sn.add_more_data(self.data_st.yt_num_list, self.data_st.yd_list)

        self.isf_func = ISFHelper()
        self.isf_func.enable_now_s_n(s_n)

    def _prepare_σμ(self, s_n, σ, μ):
        self.data_σ_μ = ΣMDataHelper(s_n, σ, μ, self.data_sn, self.data_st)
        self.data_σ_μ.hook_meta_ν_list_of_σμ_and_get_ν_st_list_inner_meta_bl()

    def _preload_self_sym_isf(self, s_n, σ, μ, ν_bl_yds):
        """这里要预读取σ，μ，ν_bl的全部分支的isf"""
        self.preload_isf_dict = {}
        for bl_ν_st in ν_bl_yds:
            _, isf_dict_info = load_isf(s_n, σ, μ, bl_ν_st, output_mode="all", ex_params=["data"])
            isf_dict = isf_dict_info["data"]
            self.preload_isf_dict[str((σ, μ, bl_ν_st))] = isf_dict

    def _get_preload_isf(self, σ, μ, ν_st):
        """因为有先前代码的检查，所以这里的get必定能取到"""
        return self.preload_isf_dict[str((σ, μ, ν_st))]


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
        # if sym_σμν == ([3, 2], [3, 1, 1], [3, 1, 1]):
        #     logger.warning("\ns_row={}, sym_key={}".format(s_row, sym_key))
        #     logger.warning("meta_σμν={}, ϵ_dict={}\n"
        #                    "meta_σμν_st={}, ϵ_st_dict={}".format(meta_σμν, ϵ_dict, meta_σμν_st, ϵ_st_dict))
        #     logger.warning("h_square={}, ϵ_4={}, ΛΛ={}".format(h_square, ϵ_4, ΛΛ))

        single_sym_isf_square = h_square * ϵ_4 * ΛΛ * single_meta_isf_square

        return single_sym_isf_square


class TestISFOverPhase(BaseCls):

    @staticmethod
    def _get_first_no_0_number_from_vector(vector):
        """提取vector中首个非0的数字和index"""
        for i, number in enumerate(vector):
            if number != 0:
                return i, number
        return None, None

    # @pytest.mark.skip("pass")
    def test_overall_phase(self):
        """
        根据书中定义，不同σμντ组合，要求ν' min下 第一个非零最小σ'μ'τ'的值为正。
        也即，不同σμντ组合，总是取第一分支的ν'，要求ντ那列第一个非零分量的值为正
        """
        for s_n in range(self.sn_min, self.sn_max + 1):
            logger.debug("s_n={}".format(s_n))
            self._prepare_sn(s_n)
            for (σ, μ) in product(self.data_sn.yd_list, repeat=2):
                _, cg_series_list = load_cg_series(s_n, σ, μ)  # 这里不使用SimpleΣMDataHelper是为了简化，不需要完整的全数据
                for ν, cg_series in zip(self.data_sn.yd_list, cg_series_list):
                    first_ν_st = self.data_sn.get_bl_yds(ν)[0]
                    for τ in range(1, cg_series + 1):
                        if cg_series == 1:
                            τ = None
                            col = ν
                        else:
                            col = (ν, τ)

                        # 到这里，合法的Sn, σ, μ, ν, ν_st_first, τ组合已经拿到
                        _, isf_dict = load_isf(s_n, σ, μ, first_ν_st,
                                               output_mode="single_col", ex_params=col)
                        isf_square_vector = isf_dict["single_col"]
                        no_0_element_index, no_0_element = self._get_first_no_0_number_from_vector(isf_square_vector)
                        # 检查
                        assert no_0_element > 0, \
                            "no_0_element={} must > 0 but not, with σ={}, μ={}, first_ν_st={}, " \
                            "col={}, no_0_element_index={}" \
                            "".format(no_0_element, σ, μ, first_ν_st, col, no_0_element_index)


class TestISF4_195b(BaseCls):

    # @pytest.mark.skip("pass")
    def test_first_isf_2_other_isf(self):
        """
        根据书中4-195b所述，利用首分支的isf可以计算出余下所有非首分支isf
        这个测试就是拿code的首分支，计算出非首分支；再拿来和code的非首分支对比
        """
        for s_n in range(self.sn_min, self.sn_max + 1):
            self._prepare_sn(s_n)
            for (σ, μ) in product(self.data_sn.yd_list, repeat=2):
                self._prepare_σμ(s_n, σ, μ)
                cg_series_list = self.data_σ_μ.cg_series_list
                for ν, cg_series in zip(self.data_sn.yd_list, cg_series_list):
                    if cg_series == 0:
                        continue
                    logger.debug("s_n={}, σ={}, μ={}, ν={}".format(s_n, σ, μ, ν))
                    # 准备
                    # 公共部分
                    ν_bl_yds = self.data_sn.get_bl_yds(ν)  # ν_st的list
                    m_ν_st_st = 1
                    # 首分支的部分(f_开头)
                    f_ν_st = ν_bl_yds[0]
                    # 到这里，首项的合法的Sn, σ, μ, ν, ν_st_first组合已经拿到
                    _, f_isf_dict_info = load_isf(s_n, σ, μ, f_ν_st, output_mode="all", ex_params=["data"])
                    f_isf_dict = f_isf_dict_info["data"]
                    f_isf_cols = f_isf_dict["cols"]
                    f_isf_rows = f_isf_dict["rows"]
                    f_isf_square_matrix = f_isf_dict["isf"]
                    # f_isf_square_list = []
                    f_isf_vector_list = []
                    for τ in range(1, cg_series + 1):
                        if cg_series == 1:
                            col = ν
                        else:
                            col = (ν, τ)
                        # 首分支的部分(f_开头)
                        f_isf_square = f_isf_square_matrix[:, f_isf_cols.index(col)]
                        f_isf_vector = Matrix([sp.sign(j) * sqrt(abs(j)) for j in f_isf_square])
                        # f_isf_square_list.append(f_isf_square)
                        f_isf_vector_list.append(f_isf_vector)

                    # 非首分支的部分(o_开头)
                    o_ν_st_list = ν_bl_yds[1::]
                    for o_ν_st in o_ν_st_list:
                        # 到这里，首项的合法的Sn, σ, μ, ν, ν_st_other组合已经拿到
                        _, o_isf_dict_info = load_isf(s_n, σ, μ, o_ν_st, output_mode="all", ex_params=["data"])
                        o_isf_dict = o_isf_dict_info["data"]
                        o_isf_cols = o_isf_dict["cols"]
                        o_isf_rows = o_isf_dict["rows"]
                        o_isf_square_matrix = o_isf_dict["isf"]
                        # o_isf_square_list = []
                        o_isf_vector_list = []
                        for τ in range(1, cg_series + 1):
                            if cg_series == 1:
                                col = ν
                            else:
                                col = (ν, τ)
                            # 首分支的部分(f_开头)
                            o_isf_square = o_isf_square_matrix[:, o_isf_cols.index(col)]
                            o_isf_vector = Matrix([sp.sign(j) * sqrt(abs(j)) for j in o_isf_square])
                            # o_isf_square_list.append(o_isf_square)
                            o_isf_vector_list.append(o_isf_vector)

                        # 下面三个m的意义分别是：
                        # 如此做的根本原因是，去掉Sn和Sn-1，f_ν和o_ν它们的ν_st_st一致，且m''一致
                        # 1 m_ν：               令当前ν_st（非ν的第一分支）的m'取1时，ν的m；
                        # 2 m_ν_st_fbl：        令当前ν_st的第一分支ν_st_st的m''取1时，ν_st_fbl（ν的第一分支）的m'
                        #                       (ν_st_st是ν_st的第一分支，可以证明它也是ν_st_fbl的分支，但不一定是第一分支)
                        # 3,m_ν_by_m_ν_st_fbl： 按照2中方法，确定m_ν_st_fbl后（也意味着确定了ν_st_fbl的杨盘），
                        #                       ν对应的m
                        ν_st_st_of_o_ν_st = self.data_st.get_bl_yds(o_ν_st)[0]  # 它也必然是f_ν_st的分支
                        m_st__o_ν_st_of_ν = self.data_st.quick_calc_m(m_ν_st_st, ν_st_st_of_o_ν_st, o_ν_st)  # 必是1
                        assert m_st__o_ν_st_of_ν == 1, "m_st__o_ν_st_of_ν={} should be 1 but not".format(
                            m_st__o_ν_st_of_ν)
                        m__o_ν = self.data_sn.quick_calc_m(m_st__o_ν_st_of_ν, o_ν_st, ν)  # 1 m_ν
                        m_st__f_ν_st_of_re = self.data_st.quick_calc_m(m_ν_st_st, ν_st_st_of_o_ν_st,
                                                                       f_ν_st)  # 2 m_ν_st_fbl
                        m__f_ν_of_f_ν_st_of_re = \
                            self.data_sn.quick_calc_m(m_st__f_ν_st_of_re, f_ν_st, ν)  # 3,m_ν_by_m_ν_st_fbl

                        # 测试
                        # 测试由首项ISF生成的非首项ISF的计算
                        # phase_vector_list = []
                        # phase_square_list = []
                        # 使用4-195b正算phase_vector
                        for f_isf_vector, o_isf_vector in zip(f_isf_vector_list, o_isf_vector_list):  # 按τ展开
                            phase_vector_element_list = []
                            for o_row in o_isf_rows:
                                flag, single_phase_vector_new = \
                                    self.isf_func._calc_isf_by_known_isf(f_isf_vector, f_isf_rows, f_ν_st,
                                                                         m_st__f_ν_st_of_re, m__f_ν_of_f_ν_st_of_re,
                                                                         o_row, o_ν_st, m_st__o_ν_st_of_ν, m__o_ν,
                                                                         ν, self.data_sn, self.data_σ_μ)
                                assert flag
                                phase_vector_element_list.append(single_phase_vector_new)
                            phase_vector = Matrix(phase_vector_element_list)
                            # phase_vector_list.append(phase_vector)
                            # phase_square = Matrix([sp.sign(i) * (i ** 2) for i in phase_vector_element_list])
                            # phase_square_list.append(phase_square)
                            assert phase_vector == o_isf_vector


class TestISFSYM(BaseCls):

    # @pytest.mark.skip("pass")
    # def test_isf_2_sym_isf(self):
    #     """
    #     根据书中4-196所述，利用对称性可以将isf对称为其他isf
    #     这个测试就是拿code的isf，计算出对称isf；再拿来和code的isf对比
    #     """
    #     pass

    @pytest.mark.skip("pass")
    def test_isf_2_self_sym_isf(self):
        """
        根据书中4-196所述，利用对称性可以将isf对称为其他isf。其中的特殊情况是对称到了自己。
        这个测试就是拿code的isf，只计算对称后还是原isf的那些；再拿来和code的isf对比
        """
        total_flag = True
        for s_n in range(self.sn_min, self.sn_max + 1):
            self._prepare_sn(s_n)
            helper = Helper(s_n - 1, s_n)
            for (σ, μ) in product(self.data_sn.yd_list, repeat=2):
                self._prepare_σμ(s_n, σ, μ)
                cg_series_list = self.data_σ_μ.cg_series_list
                for ν, cg_series in zip(self.data_sn.yd_list, cg_series_list):
                    if cg_series == 0:
                        continue
                    isf_2_self_sym_isf_flag = False
                    meta_σμν = (σ, μ, ν)
                    for current_sym_key in self.all_sym_key_list:  # 这个循环里，只要确定能对称回自己就可以了
                        sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(current_sym_key)
                        sym_σμν = tuple(meta_σμν[d] if k is False else
                                        self.data_sn.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
                        if meta_σμν == sym_σμν:
                            isf_2_self_sym_isf_flag = True
                            break
                    if isf_2_self_sym_isf_flag is False:
                        continue
                    logger.debug("\n@@@@@@@@ s_n={}, σ={}, μ={}, ν={} @@@@@@@@".format(s_n, σ, μ, ν))
                    # 到这里，确认了本组σ, μ, ν可以自对称，那么就可以先把isf读取出来，留作备用了
                    ν_bl_yds = self.data_sn.get_bl_yds(ν)  # ν_st的list
                    self._preload_self_sym_isf(s_n, σ, μ, ν_bl_yds)
                    for ν_st in ν_bl_yds:
                        meta_isf_dict = self._get_preload_isf(σ, μ, ν_st)
                        # 开始准备isf
                        meta_isf_cols = meta_isf_dict["cols"]
                        meta_isf_rows = meta_isf_dict["rows"]
                        meta_isf_square_matrix = meta_isf_dict["isf"]
                        meta_isf_square_list = []
                        # meta_isf_vector_list = []
                        for τ in range(1, cg_series + 1):
                            if cg_series == 1:
                                col = ν
                            else:
                                col = (ν, τ)
                            meta_isf_square = meta_isf_square_matrix[:, meta_isf_cols.index(col)]
                            # meta_isf_vector = Matrix([sp.sign(j) * sqrt(abs(j)) for j in meta_isf_square])
                            meta_isf_square_list.append(meta_isf_square)
                            # meta_isf_vector_list.append(meta_isf_vector)

                        # if meta_σμν == ([3, 1, 1], [3, 1, 1], [3, 1, 1]):
                        #     logger.debug("#### ν_st={}, meta_isf_square_list={}".format(ν_st, meta_isf_square_list))

                        for current_sym_key in self.all_sym_key_list:  # 这个循环里，要具体使用sym_isf了
                            sym_d3, sym_k4 = self.data_sn.get_d3_k4_by_ϵ_key(current_sym_key)
                            sym_σμν = tuple(meta_σμν[d] if k is False else
                                            self.data_sn.get_tilde(meta_σμν[d]) for d, k in zip(sym_d3, sym_k4))
                            if meta_σμν != sym_σμν:
                                continue
                            sym_isf_list = []  # (sym_ν_st, sym_row, sym_isf_square_element_list)
                            for single_meta_row in meta_isf_rows:
                                meta_row_σ, meta_row_μ, meta_row_τ = single_meta_row if len(
                                    single_meta_row) == 3 \
                                    else (*single_meta_row, None)
                                meta_σμν_st = (meta_row_σ, meta_row_μ, ν_st)
                                sym_σμν_st = tuple(meta_σμν_st[d] if k is False else
                                                   self.data_st.get_tilde(meta_σμν_st[d]) for d, k in
                                                   zip(sym_d3, sym_k4))
                                single_sym_row = (sym_σμν_st[0], sym_σμν_st[1]) if meta_row_τ is None \
                                    else (sym_σμν_st[0], sym_σμν_st[1], meta_row_τ)
                                sym_isf_dict = self._get_preload_isf(σ, μ, sym_σμν_st[2])
                                # 开始准备isf
                                sym_isf_cols = sym_isf_dict["cols"]
                                sym_isf_rows = sym_isf_dict["rows"]
                                sym_isf_square_matrix = sym_isf_dict["isf"]
                                sym_isf_square_element_list = []
                                for τ in range(1, cg_series + 1):
                                    if cg_series == 1:
                                        col = ν
                                    else:
                                        col = (ν, τ)
                                    sym_isf_square_element = \
                                        sym_isf_square_matrix[sym_isf_rows.index(single_sym_row),
                                                              sym_isf_cols.index(col)]  # 这里meta和sym是共享Sn的
                                    sym_isf_square_element_list.append(sym_isf_square_element)
                                    # sym_isf_vector_element = \
                                    #     sp.sign(sym_isf_square_element) * sp.sqrt(sym_isf_square_element)
                                sym_isf_list_element = (sym_σμν_st[2], single_sym_row, sym_isf_square_element_list)
                                sym_isf_list.append(sym_isf_list_element)

                            # if meta_σμν == ([3, 1, 1], [3, 1, 1], [3, 1, 1]) and current_sym_key == "νμσ":
                            #     logger.warning("#### sym_isf_list={}"
                            #                    "".format(sym_isf_list))

                            # 开始测试
                            σμν_ν_st_current_sym_key_flag = True
                            σμν_ν_st_current_sym_key_msg = self.isf_2_self_sym_isf_head
                            for final_row_i, final_τ in \
                                    product(range(len(sym_isf_list)), range(len(meta_isf_square_list))):
                                if cg_series == 1:
                                    τ = None  # 它也是sym_τ
                                else:
                                    τ = final_τ + 1
                                this_τ_meta_isf_square_element = meta_isf_square_list[final_τ][final_row_i]
                                meta_row = meta_isf_rows[final_row_i]
                                sym_ν_st, sym_row, sym_isf_square_list_db = sym_isf_list[final_row_i]
                                sym_isf_square_element_db = sym_isf_square_list_db[final_τ]
                                sym_isf_element = \
                                    helper.calc_sym_isf_element(sym_row, sym_σμν, τ, sym_ν_st,
                                                                current_sym_key, sym_d3, sym_k4,
                                                                self.data_sn, self.data_st,
                                                                this_τ_meta_isf_square_element)

                                # if meta_σμν == ([3, 1, 1], [3, 1, 1], [3, 1, 1]) and current_sym_key == "νμσ":
                                #     # logger.debug("$$$$ sym_isf_element={} with param={}"
                                #     #              "".format(sym_isf_element,
                                #     #                        (sym_row, sym_σμν, τ, sym_ν_st,
                                #     #                         current_sym_key, sym_d3, sym_k4,
                                #     #                         self.data_sn, self.data_st,
                                #     #                         this_τ_meta_isf_square_element)))
                                #     logger.debug("$$$$ final_row_i={}, final_τ={}, "
                                #                  "meta_row={}, sym_ν_st={}, sym_row={}, "
                                #                  "sym_isf_element={}"
                                #                  "".format(final_row_i, final_τ,
                                #                            meta_row, sym_ν_st, sym_row,
                                #                            sym_isf_element))

                                # 收集结果
                                if sym_isf_square_element_db != sym_isf_element:
                                    σμν_ν_st_current_sym_key_flag = False
                                    total_flag = False
                                    msg = self.isf_2_self_sym_isf_format.format(
                                        str(ν_st), str(meta_row), str(this_τ_meta_isf_square_element), str(τ),
                                        str(sym_ν_st), str(sym_row), str(sym_isf_element),
                                        str(sym_isf_square_element_db))
                                    σμν_ν_st_current_sym_key_msg += msg
                            if σμν_ν_st_current_sym_key_flag is False:
                                # assert sym_isf_square_element_db == sym_isf_element, \
                                #     "{}_{} ={}=> {}_{} τ={} " \
                                #     "with meta_row={}, meta_isf={}, sym_row={}, sym_isf={} " \
                                #     "and sym_isf_square_element_db={}" \
                                #     "".format(meta_σμν, ν_st, current_sym_key, sym_σμν, sym_ν_st, τ,
                                #               meta_row, this_τ_meta_isf_square_element,
                                #               sym_row, sym_isf_element,
                                #               sym_isf_square_element_db)
                                logger.warning("=={}==>".format(current_sym_key))
                                logger.info("{}".format(σμν_ν_st_current_sym_key_msg))
        assert total_flag is True, "see log"
