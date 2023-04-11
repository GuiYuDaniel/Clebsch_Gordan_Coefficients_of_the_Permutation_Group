# -*- coding:utf8 -*-
"""
验收
young_table有关的db和性质：
1，Λ^σ_mσ * Λ^σ~_mσ~ = c（常数）
"""


import pytest
from itertools import product, combinations, combinations_with_replacement, chain, permutations
from conf.cgc_config import min_s_n_of_young_table
from core.young_diagrams import load_young_diagrams, calc_young_diagram_tilde
from core.young_tableaux import load_young_table_phase_factor, get_young_tableaux_finish_s_n
from core.cg_series import load_cg_series
from utils.log import get_logger


logger = get_logger(__name__)


class TestSYM(object):

    def setup_class(self):
        # 它不改变数据库，所以不需要保护
        # 准备
        self.sn_min = max(min_s_n_of_young_table, 1)
        _, finish_s_n = get_young_tableaux_finish_s_n()
        self.sn_max = min(finish_s_n, 9)
        self.yd_list = None
        self.yd_tilde_list = None
        self.Λ_list_list = None

    def teardown_class(self):
        pass

    def _prepare_sn(self, s_n):
        _, yd_list = load_young_diagrams(s_n)
        yd_tilde_list = []
        Λ_list_list = []
        for yd in yd_list:
            _, yd_tilde = calc_young_diagram_tilde(yd, is_check_yd=False)
            yd_tilde_list.append(yd_tilde)
            _, phase_factor_list = load_young_table_phase_factor(s_n, yd, is_flag_true_if_not_s_n=False)
            Λ_list_list.append(phase_factor_list)
        self.yd_list = yd_list
        self.yd_tilde_list = yd_tilde_list
        self.Λ_list_list = Λ_list_list

    def _get_tilde(self, yd):
        return self.yd_tilde_list[self.yd_list.index(yd)]

    def _get_phase_factor_list(self, yd):
        return self.Λ_list_list[self.yd_list.index(yd)]

    def test_phase_factor_dot_self_tilde_phase_factor(self):
        """
        一个young_table的相因子，乘以它自己共轭的相因子，应该为常数，与m无关
        它的一个等价验证是：Λσ_list == reversed(Λσ~_list) or == reversed(-Λσ~_list)
        """
        for s_n in range(self.sn_min, self.sn_max + 1):
            logger.debug("s_n={}".format(s_n))
            self._prepare_sn(s_n)
            s_n_all_list = []
            # 第一个检查：Λ * Λ与m无关（都相等）
            for yd in self.yd_list:
                Λ_yd_list = self._get_phase_factor_list(yd)
                yd_tilde = self._get_tilde(yd)
                Λ_yd_tilde_list = self._get_phase_factor_list(yd_tilde)
                reverse_Λ_yd_tilde_list = list(reversed(Λ_yd_tilde_list))
                if Λ_yd_list == reverse_Λ_yd_tilde_list:
                    logger.debug("yd={}, all(Λ^σ_mσΛ^σ~_mσ~)={}, {}".format(yd, 1, "+"))  # 通过
                    s_n_all_list.append("+")
                elif Λ_yd_list == [-i for i in reverse_Λ_yd_tilde_list]:
                    logger.debug("yd={}, all(Λ^σ_mσΛ^σ~_mσ~)={}, {}".format(yd, -1, "-"))  # 通过
                    s_n_all_list.append("-")
                else:
                    # 不通过
                    assert False, "yd={}; yd_tilde={}, Λ_yd_tilde_list={}" \
                                  "Λ_yd_list={} != reverse_Λ_yd_tilde_list={} or -reverse_Λ_yd_tilde_list={}" \
                                  "".format(yd, yd_tilde, Λ_yd_tilde_list,
                                            Λ_yd_list, reverse_Λ_yd_tilde_list, [-i for i in reverse_Λ_yd_tilde_list])
            # 第二个检查：σ的Λ*Λ 与 σ~的Λ*Λ 也相等
            for yd in self.yd_list:
                rst_yd = s_n_all_list[self.yd_list.index(yd)]
                yd_tilde = self._get_tilde(yd)
                rst_yd_tilde = s_n_all_list[self.yd_list.index(yd_tilde)]
                assert rst_yd == rst_yd_tilde, "{} of yd={} must == {} of yd_tilde={}" \
                                               "".format(rst_yd, yd, rst_yd_tilde, yd_tilde)

            # 把总结果，打印一下，DEBUG等级目前就可以
            nice_msg = "\nFor Sn={}, all(Λ^σ_mσ * Λ^σ~_mσ~) is:\n".format(s_n)
            for yd, i in zip(self.yd_list, s_n_all_list):
                msg = "{:<20}: {:<5}\n".format(str(yd), i)
                nice_msg += msg
            logger.debug(nice_msg)
