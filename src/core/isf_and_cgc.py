# -*- coding: utf-8 -*-
"""
this code for creating ISF and CGC
采用对称性方法
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）本段代码主要依据置换群CGC的第二种递推计算方法
# English for 《Group Representation Theory for Physicists》 Jin-Quan Chen, Jialun Ping & Fan Wang (World Scientific)
# ISF：第四章第19节 公式：4-192（用Sn-1的本征方程计算） 4-189bc（用Sn的CG系数计算） 4-196对称性性质
# CGC：第四章第11、12节
# 整体思路是，先求Sn的ISF，再利用Sn的ISF以及Sn-1的CG，拼出Sn的CG


import copy
import json
import os
import time
# import numpy as np
import sympy as sp
from functools import lru_cache, partial
from itertools import product, combinations, combinations_with_replacement, chain, permutations
from conf.cgc_config import default_s_n, group_d3, group_k4
from conf.cgc_config import min_s_n_of_isf, min_s_n_of_cgc, min_s_n_of_branching_law, min_s_n_of_ϵ
from core.young_diagrams import load_young_diagrams, calc_young_diagram_tilde
from core.young_tableaux import load_young_table_num, load_young_table_phase_factor
from core.branching_laws import load_branching_law
from core.yamanouchi_matrix import load_yamanouchi_matrix
from core.cg_series import load_cg_series
from core.eigenvalues import load_eigenvalues
from core.symmetry_combination import load_meta_σμν, load_sym_σμν, calc_ϵ_map_dicts
from core.cgc_utils.cgc_db_typing import ISFInfo, CGCInfo, EInfo
from core.cgc_utils.cgc_local_db import get_isf_file_name, get_isf_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_cgc_file_name, get_cgc_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_ϵ_file_name, get_ϵ_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


def _debug_condition(data_σ_μ, ν_st=None, ν=None):
    """ TODO delete it """
    cond_0 = None
    cond_1 = (data_σ_μ.σ == [3, 2] and data_σ_μ.μ == [3, 1, 1] and ν_st == [2, 1, 1])
    cond_2 = (data_σ_μ.σ == [2, 1] and data_σ_μ.μ == [2, 1])
    cond_3 = (sum(data_σ_μ.σ) in [2, 3])
    cond_21 = (data_σ_μ.σ == [3, 2] and data_σ_μ.μ == [3, 1, 1] and ν_st == [3, 1])
    cond_23 = (data_σ_μ.σ == [3, 2] and data_σ_μ.μ == [3, 1, 1] and ν_st == [2, 1, 1])
    cond_27 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 1, 1] and ν_st == [3, 1])
    cond_28 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 1, 1] and ν_st == [2, 2])
    cond_29 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 1, 1] and ν_st == [2, 1, 1])
    cond_30 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [2, 2, 1] and ν_st == [3, 1])
    cond_32 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [2, 2, 1] and ν_st == [2, 1, 1])
    condition = any([cond_21, cond_23, cond_27, cond_28, cond_29, cond_30, cond_32])
    cond_4 = (data_σ_μ.σ == [2, 1, 1] and data_σ_μ.μ == [2, 1, 1] and ν_st == [1, 1, 1])
    cond_5 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 2] and ν_st == [2, 1, 1])
    cond_5_1 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 2] and ν_st in ([3, 1], [2, 1, 1]))
    cond_5_2 = (data_σ_μ.σ == [3, 2] and data_σ_μ.μ == [3, 1, 1] and ν_st in ([3, 1], [2, 1, 1]))
    cond_6 = (data_σ_μ.σ == [3, 1] and data_σ_μ.μ == [3, 1])
    cond_7 = (data_σ_μ.σ == [3, 1] and data_σ_μ.μ == [2, 1, 1])
    cond_8 = (data_σ_μ.σ == [3] and data_σ_μ.μ == [2, 1])
    cond_9 = (data_σ_μ.σ == [3, 1, 1] and data_σ_μ.μ == [3, 1, 1])
    cond_10 = (sum(data_σ_μ.σ) in [4])

    md1 = (data_σ_μ.σ == [3] and data_σ_μ.μ == [2, 1])
    md2 = cond_7 and ν == [3, 1]
    md3 = (md2 and ν == [2, 2])
    md4 = (md1 and ν == [2, 1])
    md5 = (md2 and ν_st == [2, 1])
    md6 = (data_σ_μ.σ == [4, 1] and data_σ_μ.μ == [4, 1])
    md7 = (data_σ_μ.σ == [5, 1, 1] and data_σ_μ.μ == [4, 2, 1] and ν_st == [3, 1, 1, 1])
    md8 = (data_σ_μ.σ == [3, 2, 1] and data_σ_μ.μ == [3, 2, 1])
    if md2:
        return True
    return False


def _debug_isf_cond(σ, μ, ν_st=None):
    """ TODO delete it """
    cond_0 = None
    cond_1 = (σ == [3, 1] and μ == [2, 2] and ν_st == [1, 1, 1])
    cond_2 = (σ == [3, 1] and μ == [2, 2])
    if cond_2:
        return True
    return False


def _debug_σμν_cond(σ, μ, ν=None):
    """ TODO delete it """
    cond_0 = None
    cond_1 = (σ == [3, 1] and μ == [2, 2] and ν == [2, 1, 1])
    if cond_1:
        return True
    return False


# if _debug_condition(meta_data_σ_μ):
#     logger.warning("@@@@ meta_data_σ_μ={}".format(meta_data_σ_μ))


def create_isf_and_cgc(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储isf和cgc实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    # TODO restart
    if not isinstance(s_n, int) or s_n < min(min_s_n_of_isf, min_s_n_of_cgc):
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n),
                                                                     min(min_s_n_of_isf, min_s_n_of_cgc))
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_isf_and_cgc get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, isf_finish_s_n = get_isf_finish_s_n()
    if not flag:
        err_msg = "get_isf_finish_s_n meet error with msg={}".format(isf_finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    flag, cgc_finish_s_n = get_cgc_finish_s_n()
    if not flag:
        err_msg = "get_cgc_finish_s_n meet error with msg={}".format(cgc_finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if isf_finish_s_n == min_s_n_of_isf - 1 and cgc_finish_s_n == min_s_n_of_cgc - 1:  # 特殊情况
        finish_s_n = cgc_finish_s_n
    elif isf_finish_s_n != cgc_finish_s_n:
        err_msg = "util restart code enable, isf_finish_s_n={} must eq cgc_finish_s_n={}".format(isf_finish_s_n,
                                                                                                 cgc_finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    else:
        finish_s_n = isf_finish_s_n
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} isf had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc isf s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 正式计算前的准备
    isf_func = ISFHelper()
    cgc_func = CGCHelper()
    ϵ_func = EHelper()
    ϵ_key2groups_dict, groups2ϵ_key_dict = calc_ϵ_map_dicts()

    # 开启循环计算ISF和CGC(本质上是8个参量σ μ ν τ σ’ μ’ ν’ τ’的循环，就看怎么写更优了)
    data_si = None
    # data_st = None
    # 按照从小到大的顺序，逐个计算s_i的eigenvalues并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        isf_speed_time_s_i = 0
        ϵ_speed_time_s_i = 0
        cgc_speed_time_s_i = 0
        s_i_start_time = time.time()
        logger.debug("Si={}".format(s_i))

        if s_i == 1:
            # 特殊情况，Sn=1时，没有ISF，只有CGC
            σ, μ, ν, τ, m = [1], [1], [1], None, 1
            cgc_square_dict = {(1, 1,): 1}
            cgc_speed_time = int(time.time() - s_i_start_time)
            flag, msg = save_cgc(s_i, σ, μ, ν, τ, m, cgc_square_dict, cgc_speed_time)
            if not flag:
                err_msg = "save_cgc fail by s_i={}, σ={}, μ={}, ν={}, τ={}, m={}, cgc_square_dict={} with " \
                          "msg={}".format(s_i, σ, μ, ν, τ, m, cgc_square_dict, msg)
                logger.error(err_msg)
                return False, err_msg
            s_i_cgc_speed_time = int(time.time() - s_i_start_time)
            flag, msg = save_cgc_finish_s_n(s_i, s_i_cgc_speed_time, is_check_add_one=True)
            if not flag:
                err_msg = "save_isf_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
                logger.error(err_msg)
                return False, err_msg
            continue

        isf_func.enable_now_s_n(s_i)
        cgc_func.enable_now_s_n(s_i)
        ϵ_func.enable_now_s_n(s_i)

        # 对于St，如果可以从上一次循环继承，则继承；否则要计算
        data_st = data_si if data_si is not None else DataHelper(s_i - 1, ϵ_key2groups_dict, groups2ϵ_key_dict)
        # 对于Sn=1、2，是不启动add_more_data的，后面的函数也不会调用它
        if data_si is None and s_i - 1 - 1 >= min_s_n_of_branching_law:
            data_st_st_yt_num_list = []
            _, yd_st_st_list = load_young_diagrams(s_i - 1 - 1)
            for yd_st_st in yd_st_st_list:
                _, total_num = load_young_table_num(s_i - 1 - 1, yd_st_st)
                data_st_st_yt_num_list.append(total_num)
            data_st.add_more_data(data_st_st_yt_num_list, yd_st_st_list)
        # 对于Sn，需要计算本轮的准备数据
        data_si = DataHelper(s_i, ϵ_key2groups_dict, groups2ϵ_key_dict)
        data_si.add_more_data(data_st.yt_num_list, data_st.yd_list)
        '''
        主要思路：
        a1, 计算元σμν组合首分支的ISF（ISF_matrix法）（τ>2，优先保证len4对称μσν+μ~σ~ν）
        a2, 根据情况，选择合适的计算非首分支方法（4-195c，对称，isf_matrix）
        a3, 计算元σμν的CGC
        a4, 判断和计算元σμν的ϵ（最多A2 * A2 * A3=len24，也可能比它小）
        
        b1, 计算元σμν满足的新对称组合σ*μ*ν*
        b2, 判断和计算新对称组合σ*μ*ν*的ϵ*（不同于元σμν，但可以简单分析）
        b3, 计算新对称组合σ*μ*ν*的ISF*
        b4, 计算新对称组合σ*μ*ν*的CGC*
        '''

        # 1.1，拿到Sn下所有元σμν组合
        # 1.2，过滤出元组合中所有不同的σμ（σμ_combination_of_meta_σμν_tuple_list）
        # 与σμ对应的ν列表（ν_combination_of_meta_σμν_tuple_list），要求σμ+ν必须是元组合
        # 【σμ循环开启】
        for (σ, μ) in data_si.σμ_combination_of_meta_σμν_tuple_list:
            # 初始化σμ组合数据
            data_σ_μ = ΣMDataHelper(s_i, σ, μ, data_si, data_st)
            # 2，拿到所有符合本轮ν_combination_of_meta_σμν_tuple_list的ν'
            data_σ_μ.hook_meta_ν_list_of_σμ_and_get_ν_st_list_inner_meta_bl()
            # logger.debug("now circling meta σ={}, μ={} with inner meta ν_list={}"
            #              "".format(σ, μ, data_σ_μ.meta_ν_list_of_σμ))
            logger.info("now circling meta_σμν=(σ={}, μ={}, ν_list={})".format(σ, μ, data_σ_μ.meta_ν_list_of_σμ))

            # 【第一次ν'循环开启】
            for ν_st in data_σ_μ.ν_st_list_inner_meta_bl:  # 只循环与元σμν有关的ν'
                logger.info("ν_st={}".format(ν_st))
                # 2.1，补充ISF_σμν'（注意，不是补完！）
                # 2.1.1，先尝试load一次，如果发现完成了，则跳过；否则要新建或补完
                flag, old_isf_info_dict = load_isf(s_i, σ, μ, ν_st, output_mode="all", ex_params=["data", "flags"],
                                                   is_flag_true_if_not_s_n=True)
                if not flag:
                    err_msg = "load_isf meet error with s_n={}, σ={}, μ={}, ν_st={}, " \
                              "msg={}".format(s_i, σ, μ, ν_st, old_isf_info_dict)
                    logger.error(err_msg)
                    return False, err_msg
                if old_isf_info_dict and old_isf_info_dict.get("flags", {}).get("finish_cols", None) is None:
                    continue  # 完成的标志

                # 2.1.2，计算row_index_list，如果它为[]，表示该组合下不存在ISF（理论上不应该，报错）
                isf_row_index_list = isf_func.calc_row_indexes(data_σ_μ, ν_st, data_st)
                if not isf_row_index_list:
                    err_msg = "isf_row_index_list={} must real list, pls check with data_σ_μ={}, ν_st={}, data_st={}" \
                              "".format(isf_row_index_list, data_σ_μ, ν_st, data_st)
                    logger.error(err_msg)
                    return False, err_msg

                # 2.1.3，正式补完σμν'的ISF（使用4-189a和4-195a）
                # (其中，只有元σμν是真计算；其他部分，有些被先前的元所对称了；还有一些，是本轮元的对称，需要等待对称完成，才能完整
                # （例如，[52]*[4111]=[311]是元，[52]*[4111]=[322]是它本身的对称）)
                '''
                isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                                   "cols":[[ν], ([ν], τ), ...],  # 有自由度tuple，无自由度list
                                   "isf": isf_square_matrix}  # sp.Matrix([len(rows), len(cols)])，且len(rows)=len(cols)
                '''
                meta_isf_start_time = time.time()
                flag, msg = isf_func.calc_meta_isf_dict_include_save(data_σ_μ, ν_st, isf_row_index_list,
                                                                     data_si, data_st, old_isf_info_dict)
                if not flag:
                    err_msg = "calc_meta_isf_dict_include_save meet error by data_σ_μ={}, ν_st={}, " \
                              "isf_row_index_list={}, data_si={}, data_st={}, old_isf_info_dict={} with msg={}" \
                              "".format(data_σ_μ, ν_st, isf_row_index_list, 
                                        data_si, data_st, old_isf_info_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg
                meta_isf_speed_time = int(time.time() - meta_isf_start_time)
                isf_speed_time_s_i += meta_isf_speed_time  # 总时间不需要计入old_adding_time，它在计算时已经计入了

                # 2.2，只计算元σμν部分的CGC
                # （注意，循环内，一次只能计算出ν'所对应的部分m，
                # 但是2循环结束时，元σμν对应的所有m，都应该被计算过，且必然在2中完成，不能通过对称性得到）
                meta_cgc_start_time = time.time()
                flag, msg = cgc_func.calc_cgc_by_isf_include_save(data_σ_μ, ν_st, data_si, data_st)
                if not flag:
                    err_msg = "calc_cgc_by_isf_include_save fail by data_σ_μ={}, ν_st={}, data_si={}, data_st={} " \
                              "with msg={}".format(data_σ_μ, ν_st, data_si, data_st, msg)
                    logger.error(err_msg)
                    return False, err_msg
                meta_cgc_speed_time = int(time.time() - meta_cgc_start_time)
                cgc_speed_time_s_i += meta_cgc_speed_time
            # 【第一次ν'循环结束】
            # 此时，所有σμ开头的元σμν的ISF、CGC都应该全部计算出来了。也只计算了它们。
            # (非元的：被其他元σμν轮所对称的，必然在之前被计算好了；未算的是本轮的元σμν所对称的部分，本轮下面会有计算，不在上面关于元的计算里)

            ϵ_start_time = time.time()
            meta_ν_τ_list = data_σ_μ.get_meta_ν_τ_list()
            # 【ν, τ循环开启】
            for (ν, τ) in meta_ν_τ_list:  # 这里的ν一定是与当前σμ构成元组合的ν
                # 3，只计算元σμν+τ的ϵ
                _, is_calc_ed = is_ϵ_exist(s_i, σ, μ, ν, τ)
                if is_calc_ed:
                    err_msg = "the meta combination of σ={}, μ={}, ν={}, τ={} is first appearing, " \
                              "cannot be calc ed before, pls check".format(σ, μ, ν, τ)
                    logger.error(err_msg)
                    return False, err_msg
                # 计算
                ϵ_dict, ϵ_flags = ϵ_func.calc_meta_ϵ_dict(data_σ_μ, ν, τ, data_si)
                if not ϵ_dict:
                    err_msg = "calc_meta_ϵ_dict fail by data_σ_μ={}, ν={}, τ={}, data_si={} " \
                              "with msg={}".format(data_σ_μ, ν, τ, data_si, ϵ_flags)
                    logger.error(err_msg)
                    return False, err_msg
                # 登记
                data_σ_μ.register_ϵ_dict_and_flags(ν, τ, ϵ_dict, ϵ_flags)
                # 保存
                flag, msg = save_ϵ(s_i, σ, μ, ν, τ, ϵ_dict, ϵ_flags)
                if not flag:
                    err_msg = "save_ϵ fail by s_i={}, σ={}, μ={}, ν={}, τ={}, ϵ_dict={}, ϵ_flags={} with " \
                              "msg={}".format(s_i, σ, μ, ν, τ, ϵ_dict, ϵ_flags, msg)
                    logger.error(err_msg)
                    return False, err_msg
                # 4，计算元σμν+τ的所有对称的ϵ
                flag, msg = ϵ_func.calc_symmetry_ϵ_by_meta_include_save(ϵ_dict, ϵ_flags, data_σ_μ, ν, τ, data_si)
                if not flag:
                    err_msg = "calc_symmetry_ϵ_by_meta_include_save fail by ϵ_dict={}, ϵ_flags={}, " \
                              "data_σ_μ={}, ν={}, τ={}, data_si={} " \
                              "with msg={}".format(ϵ_dict, ϵ_flags, data_σ_μ, ν, τ, data_si, msg)
                    logger.error(err_msg)
                    return False, err_msg
            ϵ_speed_time = int(time.time() - ϵ_start_time)
            ϵ_speed_time_s_i += ϵ_speed_time
            # 【ν, τ循环结束】
            # 在上述(ν, τ)循环结束时，所有σμ开头的和之前循环过的σμ，的元组合以及它们的对称组合，的ϵ都应该全部计算出来了
            # （即使有一些CGC还未计算，但它们必然不是元组合，所以只要元组合的CGC有了，ϵ也可以对称出来）

            # 5，计算元组合的对称ISF（不是整表，而是部分列）
            sym_isf_start_time = time.time()
            flag, msg = isf_func.calc_σμ_symmetry_isf_include_save(data_σ_μ, data_si, data_st)
            if not flag:
                err_msg = "calc_σμ_symmetry_isf_include_save fail by data_σ_μ={}, data_si={}, data_st={} with " \
                          "msg={}".format(data_σ_μ, data_si, data_st, msg)
                logger.error(err_msg)
                return False, err_msg
            sym_isf_speed_time = int(time.time() - sym_isf_start_time)
            isf_speed_time_s_i += sym_isf_speed_time

            # 6，计算对称性是ν结尾的CGC（σμν，σ~μ~ν，σ~μν~，σμ~ν~，μσν，μ~σ~ν，μ~σν~，μσ~ν~）
            sym_cgc_start_time = time.time()
            meta_ν_τ_list = data_σ_μ.get_meta_ν_τ_list()
            for (ν, τ) in meta_ν_τ_list:

                # if _debug_condition(data_σ_μ) and ν == [2, 1]:
                #     logger.warning("@@@@ calcing σ,μ,ν,τ={},{},{},{} and same ν sym".format(σ, μ, ν, τ))

                flag, msg = cgc_func.calc_sym_cgc_by_same_end_cgc_include_save(data_σ_μ, ν, τ, data_si)
                if not flag:
                    err_msg = "calc_sym_cgc_by_same_end_cgc_include_save fail by data_σ_μ={}, ν={}, τ={}, data_si={} " \
                              "with msg={}".format(data_σ_μ, ν, τ, data_si, msg)
                    logger.error(err_msg)
                    return False, err_msg
            # 7，计算对称性非ν结尾的CGC
            # 这里使用这种方法的主要原因是，cgc是以m_ν分开储存的。所以不如用isf再造一次σ，μ结尾的cgc，再利用它们的sym_ϵ计算同σ，μ结尾的cgc
            for meta_ν in data_σ_μ.meta_ν_list_of_σμ:
                meta_σμν = (σ, μ, meta_ν)
                new_σμν_list = [(σ, meta_ν, μ), (meta_ν, μ, σ)]  # [σνμ] & [νμσ]
                for new_σμν in new_σμν_list:
                    # 跳过条件
                    closed_loop_sym_σμν_list = []  # 如果new_σμν在这里，则其自身和同结尾yd，都已经被计算过了，可以跳过
                    for key, d3, k4 in chain(ϵ_func.σμν_0, ϵ_func.σμν_1):
                        sym_σμν = tuple(meta_σμν[d] if k is False else
                                        data_si.get_tilde(meta_σμν[d]) for d, k in zip(d3, k4))
                        if sym_σμν not in closed_loop_sym_σμν_list:
                            closed_loop_sym_σμν_list.append(sym_σμν)
                    if new_σμν in closed_loop_sym_σμν_list:
                        continue

                    # if _debug_condition(data_σ_μ):
                    #     logger.warning("@@@@ calcing new_σμν={} and same ν sym".format(new_σμν))

                    # 正式计算
                    flag, msg = cgc_func.calc_sym_μ_σ_end_cgc_include_save(new_σμν, data_si, data_st)
                    if not flag:
                        err_msg = "calc_sym_μ_σ_end_cgc_include_save fail by new_σμν={}, data_si={}, data_st={} " \
                                  "with msg={}".format(new_σμν, data_si, data_st, msg)
                        logger.error(err_msg)
                        return False, err_msg
            sym_cgc_speed_time = int(time.time() - sym_cgc_start_time)
            cgc_speed_time_s_i += sym_cgc_speed_time

        flag, msg = save_isf_finish_s_n(s_i, isf_speed_time_s_i, is_check_add_one=True)
        if not flag:
            err_msg = "save_isf_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg
        flag, msg = save_ϵ_finish_s_n(s_i, ϵ_speed_time_s_i, is_check_add_one=True)
        if not flag:
            err_msg = "save_ϵ_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg
        flag, msg = save_cgc_finish_s_n(s_i, cgc_speed_time_s_i, is_check_add_one=True)
        if not flag:
            err_msg = "save_isf_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_isf_and_cgc s_n from {} to {} done, return True, finish_s_n={}, using time={}s".format(
        finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_isf(s_n: int, σ: list, μ: list, ν_st: list, isf_square_dict: dict, speed_time: int, finish_cols=None):
    """
    自动保存或者更新isf

    ISF的落盘格式为

    <CG>/isf_info/Sn/[σ]_[μ]/[ν’].pkl
    {
    "file_name": "Sn/[σ]_[μ]/[ν’]",
    "data": isf_square_dict,
    "flags": {"speed_time": speed_time}
    }

    isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                       "cols": [[ν], ([ν], τ), ...],                   # 有自由度元组，无自由度列表, len(rows)=len(cols)
                       "isf": isf_square_matrix}                       # sp.Matrix([len(rows)]) fraction

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型； (Sn)
    [σ’][μ’]表示参与内积的两个置换群构型；[ν’]表示内积后的可能构型； (Sn-1)
    beta和beta'分别对应[ν]以及[ν’]的多重性
    isf_square_dict数值是isf的平方，符号是平方前isf系数的符号
    len(rows) = len(cols)
    speed_time表示计算用时（秒）

    例如，
    <CG>/isf_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 1].pkl
    {"rows": [([3,1],[3,1]), ([3,1],[2,1,1]), ([2,1,1],[3,1]), ([2,1,1],[2,1,1])],
     "cols": [[4,1], ([3,2],1), ([3,2],2), [3,1,1]],
     "isf": sp.Matrix([[5/12, 1/2, 1/12, 0],
                       [-1/12, 0, 5/12, 1/2],
                       [-1/12, 0, 5/12, -1/2],
                       [5/12, -1/2, 1/12, 0]])}

    因为引入非ν结尾的对称性，所以会出现只算了部分列的isf
    注意：isf_square_dict应该是更新后的字典，speed_time需要在函数外部加上旧时间
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_isf:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_isf)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν_st]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν_st, type(σ), type(μ), type(ν_st))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(isf_square_dict, dict):
        err_msg = "isf_square_dict={} with type={} must be dict".format(isf_square_dict, type(isf_square_dict))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg
    if finish_cols is not None and not isinstance(finish_cols, list):
        err_msg = "finish_cols={} with type={} must be list".format(finish_cols, type(finish_cols))
        logger.error(err_msg)
        return False, err_msg

    db_info = ISFInfo(s_n)
    _, file_name = get_isf_file_name(s_n, σ, μ, ν_st)

    if finish_cols is None:  # 它也是是否完整的标志
        table = {"file_name": file_name,
                 "data": isf_square_dict,
                 "flags": {"speed_time": speed_time}}
    else:
        table = {"file_name": file_name,
                 "data": isf_square_dict,
                 "flags": {"speed_time": speed_time,
                           "finish_cols": finish_cols}}

    _, is_calc_ed = is_isf_exist(s_n, σ, μ, ν_st)
    if is_calc_ed is False:  # 新建
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table)
        if not flag:
            return flag, msg
    else:  # 更新
        flag, msg = db_info.update_by_file_name(file_name, table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(file_name, table)
        if not flag:
            return flag, msg

    return True, None


def load_isf(s_n: int, σ: list, μ: list, ν_st: list, is_flag_true_if_not_s_n=True,
             output_mode="", ex_params=None):
    """支持两种load，分别是：

    1，只输入必要参数。返回isf矩阵及行、列index（字典）；
    2.1，output_mode='single_row'。返回ex_params所指定单行的isf及列index（字典）
    2.2，output_mode='single_col'。返回ex_params所指定单列的isf及行index（字典）
    2.3，output_mode='single_isf'。返回ex_params所指定单独行和列的isf（fraction）
    2.3，output_mode='all'。返回ex_params所指定key的字典
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_isf:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_isf)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν_st]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν_st, type(σ), type(μ), type(ν_st))
        logger.error(err_msg)
        return False, err_msg
    if output_mode == "":
        pass
    elif output_mode == "single_row":
        # 输出单行，额外输入行指标([σ'], [μ'], τ') or ([σ'], [μ'])
        if not isinstance(ex_params, tuple):
            err_msg = "ex_params={} with type={} must be tuple".format(ex_params, type(ex_params))
            logger.error(err_msg)
            return False, err_msg
    elif output_mode == "single_col":
        # 输出单列，额外输入列指标[ν] or ([ν], τ)
        if not isinstance(ex_params, (list, tuple)):
            err_msg = "ex_params={} with type={} must be tuple or list".format(ex_params, type(ex_params))
            logger.error(err_msg)
            return False, err_msg
    elif output_mode == "single_isf":
        # 输出单独isf，额外输入行列指标([σ'], [μ'], τ') or ([σ'], [μ']) / [ν] or ([ν], τ)
        if not isinstance(ex_params[0], tuple):
            err_msg = "ex_params[0]={} with type={} must be tuple".format(ex_params[0], type(ex_params[0]))
            logger.error(err_msg)
            return False, err_msg
        if not isinstance(ex_params[1], (list, tuple)):
            err_msg = "ex_params[1]={} with type={} must be tuple".format(ex_params[1], type(ex_params[1]))
            logger.error(err_msg)
            return False, err_msg
    elif output_mode == "all":
        if not isinstance(ex_params, (list, tuple)):
            err_msg = "ex_params={} with type={} must be list or tuple".format(ex_params, type(ex_params))
            logger.error(err_msg)
            return False, err_msg
        if any((i not in ["data", "flags"])for i in ex_params):
            err_msg = "ex_params={} of output_mode=all must in all ['data', 'flags']".format(ex_params)
            logger.error(err_msg)
            return False, err_msg
    else:
        err_msg = "output_mode={} must in ['', 'single_row', 'single_col', 'single_isf', 'all'] but not".format(
            output_mode)
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_isf_file_name(s_n, σ, μ, ν_st)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, σ={}, μ={}, ν_st={} because {}".format(
            s_n, σ, μ, ν_st, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = ISFInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query isf with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        if output_mode == "all":
            rst_dict = {}
            for key in ex_params:
                if key in data:
                    rst_dict[key] = data[key]
            return True, rst_dict  # bingo(1/5)！
        isf_square_dict = data.get("data")
        if isinstance(isf_square_dict, dict) and isf_square_dict:
            if output_mode == "single_row":
                # 输出单行，额外输入行指标([σ'], [μ'], τ') or ([σ'], [μ'])
                ex_row = ex_params[: -1] if None in ex_params else ex_params  # 使它也能接受([σ'], [μ'], None)
                if ex_row in isf_square_dict.get("rows", []):
                    row_index = isf_square_dict.get("rows").index(ex_row)
                    rst_dict = {"cols": isf_square_dict["cols"],
                                "single_row": isf_square_dict["isf"][row_index, :]}
                    return True, rst_dict  # bingo(2/5)！
                else:
                    err_msg = "ex_row={} should in rows={} with isf_square_dict={} but not, pls check".format(
                        ex_row, isf_square_dict.get("rows"), isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            elif output_mode == "single_col":
                # 输出单列，额外输入列指标[ν] or ([ν], τ)
                ex_col = ex_params[0] if None in ex_params else ex_params  # 使它也能接受([ν], None)
                if ex_col in isf_square_dict.get("cols", []):
                    cols_index = isf_square_dict.get("cols").index(ex_col)
                    rst_dict = {"rows": isf_square_dict["rows"],
                                "single_col": isf_square_dict["isf"][:, cols_index]}
                    return True, rst_dict  # bingo(3/5)！
                else:
                    err_msg = "ex_col={} should in cols={} with isf_square_dict={} but not, pls check".format(
                        ex_col, isf_square_dict.get("cols"), isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            elif output_mode == "single_isf":
                ex_row = ex_params[0][: -1] if None in ex_params[0] else ex_params[0]  # 使它也能接受([σ'], [μ'], None)
                ex_col = ex_params[1][0] if None in ex_params[1] else ex_params[1]  # 使它也能接受([ν], None)
                if ex_row in isf_square_dict.get("rows", []) and ex_params[1] in isf_square_dict.get("cols", []):
                    row_index = isf_square_dict.get("rows").index(ex_row)
                    cols_index = isf_square_dict.get("cols").index(ex_col)
                    rst = isf_square_dict["isf"][row_index, cols_index]
                    return True, rst  # bingo(4/5)！
                else:
                    err_msg = "ex_row={} should in rows={} and ex_col={} should in cols={} " \
                              "with isf_square_dict={} but not, pls check".format(ex_row, isf_square_dict.get("rows"),
                                                                                  ex_col, isf_square_dict.get("cols"),
                                                                                  isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            else:  # output_mode = ""
                return True, isf_square_dict  # bingo(5/5)！

        else:
            err_msg = "isf_square_dict queried from db, but cannot get it from data" \
                      "with data={}, isf_square_dict={} from db".format(data, isf_square_dict)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist isf_square_dict db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def is_isf_exist(s_n: int, σ: list, μ: list, ν_st: list):
    """只关注存在性，不真正读取数据"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_isf:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_isf)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν_st]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν_st, type(σ), type(μ), type(ν_st))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_isf_file_name(s_n, σ, μ, ν_st)
    # flag, file_name = get_isf_file_name(s_n, σ, μ, ν_st, is_full_path=True)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, σ={}, μ={}, ν_st={} because {}".format(
            s_n, σ, μ, ν_st, file_name)
        logger.error(err_msg)
        return False, err_msg
    return ISFInfo(s_n).exist_by_file_name(file_name)


def save_isf_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_isf:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_isf)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_isf_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = ISFInfo(0)
    _, finish_file_name = get_isf_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {"S{}".format(s_n): s_n_speed_time}}}

    if finish_s_n_before == min_s_n_of_isf - 1:
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table, point_key="flags")
        if not flag:
            return flag, msg
    else:
        flag, data = db_info.query_by_file_name(finish_file_name)
        if not flag:
            return flag, data
        history_times = data.get("flags", {}).get("history_times")
        if not history_times or not isinstance(history_times, dict):  # 有，就不应该是空字典
            err_msg = "old history_times={} must real dict, but not, with data={}".format(history_times, data)
            logger.error(err_msg)
            return False, err_msg
        # 更新table里的历史时间
        history_times.update({"S{}".format(s_n): s_n_speed_time})
        table["flags"]["history_times"] = history_times
        flag, msg = db_info.update_by_file_name(finish_file_name, partial_table=table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(finish_file_name, partial_table=table, point_key="flags")
        if not flag:
            return flag, msg

    return True, None


def get_isf_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 0
    """
    _, finish_file_name = get_isf_finish_s_n_name()
    flag, data = ISFInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get_isf_finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 1")
        return True, min_s_n_of_isf - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_isf:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_isf, data)
        return False, err_msg


def save_cgc(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), m: int, cgc_square_data, speed_time: int=0):
    """
    这个db用来存CGC

    <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_τ_m.pkl  # 无多重性，则为[ν]_m.pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_τ_m,  # 无多重性，则为[ν]_m
    "data": cgc_square_dict,
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    τ对应[ν]的多重性;
    m是[ν]的yamanouchi序;
    cgc_square_dict数值是cgc的平方，符号是平方前cgc的符号;
    speed_time表示计算用时（秒）

    例如，
    <CG>/cgc_info/S4/[2, 2]_[3, 1]/[2, 1, 1]_m2.pkl
    {(2, 1): 1/2, (1, 3): 1/4, (2, 2): 1/4}
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν, type(σ), type(μ), type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and not isinstance(τ, int):
        err_msg = "τ={} with type={} must be None or int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(m, int) or s_n <= 0:
        err_msg = "m={} with type={} must be int and > 0".format(m, type(m))
        logger.error(err_msg)
        return False, err_msg
    # if not isinstance(cgc_square_dict, dict):
    #     err_msg = "cgc_square_dict={} with type={} must be dict".format(cgc_square_dict, type(cgc_square_dict))
    #     logger.error(err_msg)
    #     return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = CGCInfo(s_n)
    _, file_name = get_cgc_file_name(s_n, σ, μ, ν, τ, m)

    table = {"file_name": file_name,
             "data": cgc_square_data,
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def load_cgc(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), m: int, is_flag_true_if_not_s_n=True):
    """
    取得s_n下指定[σ][μ][ν] τ m 的cg系数字典
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_cgc:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_cgc)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(σ, list) or not isinstance(μ, list) or not isinstance(ν, list):
        err_msg = "type(σ={})={} type(μ={})={} type(ν={})={} must be list".format(
            σ, type(σ), μ, type(μ), ν, type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and (not isinstance(τ, int) or τ <= 0):
        err_msg = "type(τ={})={} must be None or real int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(m, int) or m <= 0:
        err_msg = "type(m={})={} must be real int".format(m, type(m))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, τ, m)
    if not flag:
        err_msg = "cannot get file_name by s_n={} σ={}, μ={}, ν={}, τ={}, m={}, because {}".format(
            s_n, σ, μ, ν, τ, m, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = CGCInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query cgc with s_n={} σ={}, μ={}, ν={}, τ={}, m={}, because {}".format(
            s_n, σ, μ, ν, τ, m, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        cgc_dict = data.get("data")
        if isinstance(cgc_dict, dict) and cgc_dict:
            # 只检查有没有 不对内容做检查了
            return True, cgc_dict  # bingo！

        else:
            err_msg = "cgc_dict queried from db, but cannot get it from data" \
                      "with data={}, cgc_dict={} from db".format(data, cgc_dict)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist cgc_dict db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


load_cgc_with_m1 = partial(load_cgc, m=1, is_flag_true_if_not_s_n=False)


def is_cgc_exist(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), m: int):
    """只关注存在性，不真正读取数据"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_isf:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_isf)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
        err_msg = "all [σ={}, μ={}, ν={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν, type(σ), type(μ), type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and (not isinstance(τ, int) or τ <= 0):
        err_msg = "type(τ={})={} must be None or real int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(m, int) or m <= 0:
        err_msg = "type(m={})={} must be real int".format(m, type(m))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, τ, m)
    # flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, τ, m, is_full_path=True)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, σ={}, μ={}, ν={}, τ={}, m={} because {}".format(
            s_n, σ, μ, ν, τ, m, file_name)
        logger.error(err_msg)
        return False, err_msg
    return CGCInfo(s_n).exist_by_file_name(file_name)


def save_cgc_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_cgc_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = CGCInfo(0)
    _, finish_file_name = get_cgc_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {"S{}".format(s_n): s_n_speed_time}}}

    if finish_s_n_before == min_s_n_of_cgc - 1:
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table, point_key="flags")
        if not flag:
            return flag, msg
    else:
        flag, data = db_info.query_by_file_name(finish_file_name)
        if not flag:
            return flag, data
        history_times = data.get("flags", {}).get("history_times")
        if not history_times or not isinstance(history_times, dict):  # 有，就不应该是空字典
            err_msg = "old history_times={} must real dict, but not, with data={}".format(history_times, data)
            logger.error(err_msg)
            return False, err_msg
        # 更新table里的历史时间
        history_times.update({"S{}".format(s_n): s_n_speed_time})
        table["flags"]["history_times"] = history_times
        flag, msg = db_info.update_by_file_name(finish_file_name, partial_table=table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(finish_file_name, partial_table=table, point_key="flags")
        if not flag:
            return flag, msg

    return True, None


def get_cgc_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 0
    """
    _, finish_file_name = get_cgc_finish_s_n_name()
    flag, data = CGCInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get_cgc_finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_cgc - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_cgc:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_cgc, data)
        return False, err_msg


def save_ϵ(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), ϵ_dict: dict, flag_dict: dict):
    """
    这个db用来存ϵ

    <CG>/ϵ_info/Sn/[σ]_[μ]/[ν]_τ.pkl  # 无多重性，则为[ν].pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_τ,  # 无多重性，则为[ν]
    "data": ϵ_dict,
    "flags": flag_dict
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    τ对应[ν]的多重性;

    例如，
    <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 2]_2.pkl
    {"data": {"σμν": 1, "μσν": 1, "σ~μ~ν": -1, "μ~σ~ν": -1,
              "σ~μν~": 1, "μ~σν~": -1, "σμ~ν~": -1, "μσ~ν~": 1},
     "flags": {"σμν": (1, 1), "μσν": (1, 1), "σ~μ~ν": (6, 6), "μ~σ~ν": (6, 6),
               "σ~μν~": (6, 4), "μ~σν~": (3, 1), "σμ~ν~": (1, 3), "μσ~ν~": (4, 6)}}
    <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[2, 1, 1]_2.pkl
    {"data": {"σμν": 1, "μσν": -1, "σ~μ~ν": -1, "μ~σ~ν": 1,
              "σ~μν~": -1, "μ~σν~": -1, "σμ~ν~": 1, "μ~σ~ν": 1},
     "flags": {"σμν": (1, 4), "μσν": (4, 1), "σ~μ~ν": (6, 3), "μ~σ~ν": (3, 6),
               "σ~μν~": (6, 1), "μ~σν~": (6, 1), "σμ~ν~": (1, 6), "μ~σ~ν": (1, 6)}}
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν, type(σ), type(μ), type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and not isinstance(τ, int):
        err_msg = "τ={} with type={} must be None or int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(ϵ_dict, dict):
        err_msg = "ϵ_dict={} with type={} must be dict".format(ϵ_dict, type(ϵ_dict))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(flag_dict, dict):
        err_msg = "flag_dict={} with type={} must be dict".format(flag_dict, type(flag_dict))
        logger.error(err_msg)
        return False, err_msg

    db_info = EInfo(s_n)
    _, file_name = get_ϵ_file_name(s_n, σ, μ, ν, τ)

    table = {"file_name": file_name,
             "data": ϵ_dict,
             "flags": flag_dict}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


# def update_ϵ(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), ϵ_dict: dict):
#     if not isinstance(s_n, int) or s_n <= 0:
#         err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
#         logger.error(err_msg)
#         return False, err_msg
#     if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
#         err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
#             σ, μ, ν, type(σ), type(μ), type(ν))
#         logger.error(err_msg)
#         return False, err_msg
#     if τ is not None and not isinstance(τ, int):
#         err_msg = "τ={} with type={} must be None or int".format(τ, type(τ))
#         logger.error(err_msg)
#         return False, err_msg
#     if not isinstance(ϵ_dict, dict):
#         err_msg = "ϵ_dict={} with type={} must be dict".format(ϵ_dict, type(ϵ_dict))
#         logger.error(err_msg)
#         return False, err_msg
#
#     db_info = EInfo(s_n)
#     _, file_name = get_ϵ_file_name(s_n, σ, μ, ν, τ)
#
#     new_table = {"file_name": file_name,
#                  "data": ϵ_dict,
#                  "flags": {}}
#     flag, msg = db_info.update_by_file_name(file_name, new_table)
#     if not flag:
#         return flag, msg
#     flag, msg = db_info.update_txt_by_file_name(file_name, new_table)
#     if not flag:
#         return flag, msg
#
#     return True, None


def load_ϵ(s_n: int, σ: list, μ: list, ν: list, τ: (int, None), is_with_flags=False, is_flag_true_if_not_s_n=True):
    """
    取得s_n下指定[σ][μ][ν] τ 的ϵ字典
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int):
        err_msg = "s_n={} with type={} must be int".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(σ, list) or not isinstance(μ, list) or not isinstance(ν, list):
        err_msg = "type(σ={})={} type(μ={})={} type(ν={})={} must be list".format(
            σ, type(σ), μ, type(μ), ν, type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and (not isinstance(τ, int) or τ <= 0):
        err_msg = "type(τ={})={} must be None or real int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_ϵ_file_name(s_n, σ, μ, ν, τ)
    if not flag:
        err_msg = "cannot get file_name by s_n={} σ={}, μ={}, ν={}, τ={}, because {}".format(
            s_n, σ, μ, ν, τ, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = EInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query ϵ with s_n={} σ={}, μ={}, ν={}, τ={}, because {}".format(
            s_n, σ, μ, ν, τ, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        ϵ_dict = data.get("data")
        if is_with_flags:
            flag_dict = data.get("flags")
            if isinstance(ϵ_dict, dict) and ϵ_dict and isinstance(flag_dict, dict) and flag_dict:
                # 只检查有没有 不对内容做检查了
                return True, {"data": ϵ_dict, "flags": flag_dict}  # bingo！

            else:
                err_msg = "ϵ_dict and flag_dict queried from db, but cannot get it from data" \
                          "with data={}, ϵ_dict={} flag_dict={} from db".format(data, ϵ_dict, flag_dict)
                logger.error(err_msg)
                return False, err_msg
        else:
            if isinstance(ϵ_dict, dict) and ϵ_dict:
                # 只检查有没有 不对内容做检查了
                return True, ϵ_dict  # bingo！

            else:
                err_msg = "ϵ_dict queried from db, but cannot get it from data" \
                          "with data={}, ϵ_dict={} from db".format(data, ϵ_dict)
                logger.error(err_msg)
                return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist ϵ_dict db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def is_ϵ_exist(s_n: int, σ: list, μ: list, ν: list, τ: (int, None)):
    """只关注存在性，不真正读取数据"""
    if not isinstance(s_n, int):
        err_msg = "s_n={} with type={} must be int".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(σ, list) or not isinstance(μ, list) or not isinstance(ν, list):
        err_msg = "type(σ={})={} type(μ={})={} type(ν={})={} must be list".format(
            σ, type(σ), μ, type(μ), ν, type(ν))
        logger.error(err_msg)
        return False, err_msg
    if τ is not None and (not isinstance(τ, int) or τ <= 0):
        err_msg = "type(τ={})={} must be None or real int".format(τ, type(τ))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_ϵ_file_name(s_n, σ, μ, ν, τ)
    if not flag:
        err_msg = "cannot get file_name by s_n={} σ={}, μ={}, ν={}, τ={}, because {}".format(
            s_n, σ, μ, ν, τ, file_name)
        logger.error(err_msg)
        return False, err_msg
    return EInfo(s_n).exist_by_file_name(file_name)


def save_ϵ_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_ϵ_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = EInfo(0)
    _, finish_file_name = get_ϵ_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {"S{}".format(s_n): s_n_speed_time}}}

    if finish_s_n_before == min_s_n_of_ϵ - 1:
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table, point_key="flags")
        if not flag:
            return flag, msg
    else:
        flag, data = db_info.query_by_file_name(finish_file_name)
        if not flag:
            return flag, data
        history_times = data.get("flags", {}).get("history_times")
        if not history_times or not isinstance(history_times, dict):  # 有，就不应该是空字典
            err_msg = "old history_times={} must real dict, but not, with data={}".format(history_times, data)
            logger.error(err_msg)
            return False, err_msg
        # 更新table里的历史时间
        history_times.update({"S{}".format(s_n): s_n_speed_time})
        table["flags"]["history_times"] = history_times
        flag, msg = db_info.update_by_file_name(finish_file_name, partial_table=table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(finish_file_name, partial_table=table, point_key="flags")
        if not flag:
            return flag, msg

    return True, None


def get_ϵ_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 0
    """
    _, finish_file_name = get_ϵ_finish_s_n_name()
    flag, data = EInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get_ϵ_finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_ϵ - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_ϵ:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_ϵ, data)
        return False, err_msg


class SimpleΣMDataHelper(object):
    """简单版的ΣMDataHelper，比之少读取一些数据"""

    def __init__(self, s_n, σ, μ, data_sn_cls, data_st_cls):
        if not isinstance(s_n, int) or s_n < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_n={} must be int and >= {}".format(s_n, min(min_s_n_of_isf, min_s_n_of_cgc)))
        self.s_n = s_n
        self.s_t = s_n - 1
        self.σ = σ
        self.μ = μ
        self.data_sn_cls = data_sn_cls
        self.data_st_cls = data_st_cls

        self.cg_series_list = None
        self.bl_yds_of_σ = None
        self.bl_yds_of_μ = None
        self.cg_series_st_list_dict = None
        self._init_data()

    def __str__(self):
        return "σ={}, μ={} data simple class".format(self.σ, self.μ)

    def _init_data(self):
        # 1, Sn、σ、μ对应的cg_series_list（yd序）
        flag, cg_series_list = load_cg_series(self.s_n, self.σ, self.μ, is_flag_true_if_not_s_n=False)  # [0, 1, 0]
        if not flag:
            err_msg = "get cg_series_list meet error with s_i={}, σ={}, μ={}, " \
                      "msg={}".format(self.s_n, self.σ, self.μ, cg_series_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.cg_series_list = cg_series_list

        # 2, σ、μ对应的分支律合法yds
        self.bl_yds_of_σ = self.data_sn_cls.get_bl_yds(self.σ)
        self.bl_yds_of_μ = self.data_sn_cls.get_bl_yds(self.μ)

        # 3, （tuple(σ'), tuple(μ')）为key， cg_series_st_list为value的dict
        flag, cg_series_st_list_dict = self._load_cg_series_list_dict_by_combination_of_bl_yds()
        if not flag:
            err_msg = "get cg_series_st_list_dict meet error with msg={}".format(cg_series_st_list_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.cg_series_st_list_dict = cg_series_st_list_dict

    def _load_cg_series_list_dict_by_combination_of_bl_yds(self):
        """根据σ和μ在分支律下的有限种[σ’] [μ’]组合，一次性取出所有cg_series，按照dict形式返回"""
        yd_st_list = self.data_st_cls.yd_list
        cg_series_st_list_dict = {}
        for bl_σ, bl_μ in product(self.bl_yds_of_σ, self.bl_yds_of_μ):
            flag, cg_series_st_list = load_cg_series(self.s_t, bl_σ, bl_μ, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cg_series_st_list for _load_cg_series_list_dict_by_combination_of_bl_yds " \
                          "meet error with s_t={}, bl_σ={}, bl_μ={}, " \
                          "msg={}".format(self.s_t, bl_σ, bl_μ, cg_series_st_list)
                logger.error(err_msg)
                return False, err_msg

            if len(cg_series_st_list) != len(yd_st_list):
                err_msg = "len={} of cg_series_st_list={} must eq len={} of yd_st_list={}".format(
                    len(cg_series_st_list), cg_series_st_list, len(yd_st_list), yd_st_list)
                logger.error(err_msg)
                return False, err_msg
            cg_series_st_list_dict[(tuple(bl_σ), tuple(bl_μ))] = cg_series_st_list
        return True, cg_series_st_list_dict

    def get_cg_series_st_list(self, σ_st, μ_st):
        key = (tuple(σ_st), tuple(μ_st))
        return self.cg_series_st_list_dict.get(key)

    def get_cg_series_st(self, σ_st, μ_st, ν_st, yd_st_list):
        cg_series_st_list = self.get_cg_series_st_list(σ_st, μ_st)
        return cg_series_st_list[yd_st_list.index(ν_st)]


class ΣMDataHelper(SimpleΣMDataHelper):
    """这里封装供模块内部，σ、μ循环应该取得的数据"""

    def __init__(self, s_n, σ, μ, data_sn_cls, data_st_cls):
        if not isinstance(s_n, int) or s_n < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_n={} must be int and >= {}".format(s_n, min(min_s_n_of_isf, min_s_n_of_cgc)))
        super(ΣMDataHelper, self).__init__(s_n, σ, μ, data_sn_cls, data_st_cls)
        self.in_matrix_σ_dict = None
        self.in_matrix_μ_dict = None
        self._add_init_data()

        self.meta_ν_list_of_σμ = []
        self.ν_st_list_inner_meta_bl = []
        self.real_calc_isf_ν_st_list = []
        self.isf_square_dict_dict_by_ν_st = {}
        self.isf_finish_cols_dict_by_ν_st = {}
        self.isf_meta_ν_τ_list = []
        self.ϵ_dict_dict = {}
        self.ϵ_flags_dict = {}

    def __str__(self):
        return "σ={}, μ={} data class".format(self.σ, self.μ)

    def _add_init_data(self):
        # TODO, 1, 这个循环也许可以省略掉了（当彻底确立4-193c的时候）
        # TODO, 2, 那时，这个load也可以放在DataSn里了
        # 4, (in) matrix，计算的越早，重复的次数越少
        in_matrix_σ_dict = {}
        for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
            in_key = (i, self.s_n,)
            flag, in_matrix_σ = load_yamanouchi_matrix(self.s_n, self.σ, in_key, mode="in",
                                                       is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get in_matrix_σ with s_n={}, σ={}, in_key={} meet error with " \
                          "msg={}".format(self.s_n, self.σ, in_key, in_matrix_σ)
                logger.error(err_msg)
                raise Exception(err_msg)
            in_matrix_σ_dict[in_key] = in_matrix_σ
        self.in_matrix_σ_dict = in_matrix_σ_dict

        if self.μ == self.σ:
            self.in_matrix_μ_dict = in_matrix_σ_dict
        else:
            in_matrix_μ_dict = {}
            for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
                in_key = (i, self.s_n,)
                flag, in_matrix_μ = load_yamanouchi_matrix(self.s_n, self.μ, in_key, mode="in",
                                                           is_flag_true_if_not_s_n=False)
                if not flag:
                    err_msg = "get in_matrix_μ with s_n={}, μ={}, in_key={} meet error with " \
                              "msg={}".format(self.s_n, self.μ, in_key, in_matrix_μ)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                in_matrix_μ_dict[in_key] = in_matrix_μ
            self.in_matrix_μ_dict = in_matrix_μ_dict
            
    def hook_meta_ν_list_of_σμ_and_get_ν_st_list_inner_meta_bl(self):
        # 挂载meta_ν_list_of_σμ
        # 得到ν_st_list_inner_meta_bl
        self.meta_ν_list_of_σμ = self.data_sn_cls.get_ν_combination_of_meta_σμν_tuple_list_by_σμ(self.σ, self.μ)
        ν_st_list_inner_meta_bl = []
        tmp_repeatable_meta_ν_bl_list = []
        for meta_ν in self.meta_ν_list_of_σμ:
            tmp_repeatable_meta_ν_bl_list += self.data_sn_cls.get_bl_yds(meta_ν)
        for ν_st in self.data_st_cls.yd_list:
            if ν_st in tmp_repeatable_meta_ν_bl_list:
                ν_st_list_inner_meta_bl.append(ν_st)
        self.ν_st_list_inner_meta_bl = ν_st_list_inner_meta_bl

    def register_real_calc_isf_ν_st_list(self, ν_st):
        self.real_calc_isf_ν_st_list.append(ν_st)

    def get_real_calc_isf_ν_st_list(self):
        return self.real_calc_isf_ν_st_list

    def register_isf_square_dict_and_finish_cols_by_ν_st(self, ν_st, isf_square_dict, finish_cols):
        self.isf_square_dict_dict_by_ν_st[tuple(ν_st)] = isf_square_dict
        self.isf_finish_cols_dict_by_ν_st[tuple(ν_st)] = finish_cols

    def get_isf_square_dict_by_ν_st(self, ν_st):
        return self.isf_square_dict_dict_by_ν_st[tuple(ν_st)]

    def get_isf_finish_cols_by_ν_st(self, ν_st):
        return self.isf_finish_cols_dict_by_ν_st[tuple(ν_st)]

    # def register_isf_ν_τ_list(self, ν, τ):
    #     if (ν, τ) not in self.isf_ν_τ_list:
    #         self.isf_ν_τ_list.append((ν, τ,))

    def register_meta_ν_τ_list(self, ν, τ):
        if (ν, τ) not in self.isf_meta_ν_τ_list and ν in self.meta_ν_list_of_σμ:
            self.isf_meta_ν_τ_list.append((ν, τ,))

    def get_meta_ν_τ_list(self):
        return self.isf_meta_ν_τ_list

    def register_ϵ_dict_and_flags(self, ν, τ, ϵ_dict, ϵ_flags):
        key = tuple(ν) if τ is None else (tuple(ν), τ)
        self.ϵ_dict_dict[key] = ϵ_dict
        self.ϵ_flags_dict[key] = ϵ_flags

    def get_ϵ_dict_and_flags(self, ν, τ):
        key = tuple(ν) if τ is None else (tuple(ν), τ)
        return self.ϵ_dict_dict.get(key, {}), self.ϵ_flags_dict.get(key, {})


class DataHelper(object):
    """这里封装供模块内部，Si循环应该取得的数据"""

    def __init__(self, s_k, ϵ_key2groups_dict, groups2ϵ_key_dict):
        if not isinstance(s_k, int) or s_k < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_k={} must be int and >= 1".format(s_k))
        self.s_n = s_k  # 既可以实例化Sn，也可以实例化St

        self.yd_list = None
        self.bl_yd_list_dict = None
        self.bl_row_list_dict = None
        self.yt_num_list = None
        self.eigenvalue_list = None
        self.yd_tilde_list = None
        self.phase_factor_list_list = None
        self.butler_phase_factor_dict_dict = None  # (英文版4-196')
        self.meta_σμν_tuple_list = None
        self.σμ_combination_of_meta_σμν_tuple_list = None
        self.ν_combination_of_meta_σμν_tuple_list_list = None
        self.ϵ_key2groups_dict = ϵ_key2groups_dict
        self.groups2ϵ_key_dict = groups2ϵ_key_dict
        self._init_data()

        # 需手动加载
        self.offset_dict_list = None
        self.yt_st_num_list = None
        self.yd_st_list = None

    def __str__(self):
        return "s_n={} data class".format(self.s_n)

    def _init_data(self):
        # 1, Sn下的yd列表
        flag, yd_list = load_young_diagrams(self.s_n, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get yd_list meet error with s_n={}, msg={}".format(self.s_n, yd_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.yd_list = yd_list

        # 2, 分支律，这里先只拿before_YD，其他根据情况可添加
        # 注意：before_YD是Sn的合法子yds，本身属于St的yd_list的子集
        bl_yd_list_dict = {}
        bl_row_list_dict = {}
        for yd in yd_list:
            flag, bl_dict = load_branching_law(self.s_n, yd, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get bl_dict meet error with s_n={}, yd={}, msg={}".format(self.s_n, yd, bl_dict)
                logger.error(err_msg)
                raise Exception(err_msg)
            bl_yd_list_dict[tuple(yd)] = bl_dict.get("before_YD")
            bl_row_list_dict[tuple(yd)] = bl_dict.get("rows")
        self.bl_yd_list_dict = bl_yd_list_dict
        self.bl_row_list_dict = bl_row_list_dict

        # 3, yd对应的yt_num
        yt_num_list = []
        for yd in yd_list:
            flag, total_num = load_young_table_num(self.s_n, yd, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get yd total_num meet error with s_i={}, yd={}, msg={}".format(self.s_n, yd, total_num)
                logger.error(err_msg)
                raise Exception(err_msg)
            yt_num_list.append(total_num)
        self.yt_num_list = yt_num_list

        # 4, eigenvalue，按照yd序
        flag, eigenvalue_list = load_eigenvalues(self.s_n, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get eigenvalue_list meet error with s_n={}, msg={}".format(self.s_n, eigenvalue_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.eigenvalue_list = eigenvalue_list

        # 5, yd对应的对称yd
        yd_tilde_list = []
        for yd in yd_list:
            flag, yd_tilde = calc_young_diagram_tilde(yd, is_check_yd=False)
            if not flag:
                err_msg = "calc yd_tilde meet error with s_n={}, yd={}, msg={}".format(self.s_n, yd, yd_tilde)
                logger.error(err_msg)
                raise Exception(err_msg)
            yd_tilde_list.append(yd_tilde)
        self.yd_tilde_list = yd_tilde_list

        # 6, yd对应的phase_factor_list
        phase_factor_list_list = []
        for yd in yd_list:
            flag, phase_factor_list = load_young_table_phase_factor(self.s_n, yd, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_young_table_phase_factor meet error with s_n={}, yd={}, msg={}".format(
                    self.s_n, yd, phase_factor_list)
                logger.error(err_msg)
                raise Exception(err_msg)
            phase_factor_list_list.append(phase_factor_list)
        self.phase_factor_list_list = phase_factor_list_list

        # 7，计算factorization property of the phase factor（Butler 1979）
        butler_phase_factor_dict_dict = {}
        for yd in self.yd_list:
            butler_phase_factor_dict = {}
            for yd_st, row_st in zip(self.bl_yd_list_dict[tuple(yd)], self.bl_row_list_dict[tuple(yd)]):
                nb = sum(yd[row_st + 1:])  # where nb is the number of boxes below the box labelled with n in yt of Sn
                if nb % 2 == 0:
                    butler_phase_factor_dict[tuple(yd_st)] = 1
                else:
                    butler_phase_factor_dict[tuple(yd_st)] = -1
            butler_phase_factor_dict_dict[tuple(yd)] = butler_phase_factor_dict
        self.butler_phase_factor_dict_dict = butler_phase_factor_dict_dict

        # 8, 提取meta_σμν_tuple_list
        flag, meta_σμν_tuple_list = load_meta_σμν(self.s_n, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "load_meta_σμν meet error with s_n={}, msg={}".format(self.s_n, meta_σμν_tuple_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.meta_σμν_tuple_list = meta_σμν_tuple_list

        # 9, 计算σμ_combination_of_meta_σμν_tuple_list和ν_combination_of_meta_σμν_tuple_list_list
        σμ_combination_of_meta_σμν_tuple_list = []
        ν_combination_of_meta_σμν_tuple_list_list = []
        for meta_σμν in meta_σμν_tuple_list:
            if meta_σμν[0:2] not in σμ_combination_of_meta_σμν_tuple_list:
                σμ_combination_of_meta_σμν_tuple_list.append(meta_σμν[0:2])
                ν_combination_of_meta_σμν_tuple_list_list.append([meta_σμν[-1]])
            else:
                ν_combination_of_meta_σμν_tuple_list_list[-1].append(meta_σμν[-1])
        self.σμ_combination_of_meta_σμν_tuple_list = σμ_combination_of_meta_σμν_tuple_list
        self.ν_combination_of_meta_σμν_tuple_list_list = ν_combination_of_meta_σμν_tuple_list_list

    def add_more_data(self, yt_st_num_list, yd_st_list):
        # 无法仅通过Sn，还需要传入其他参数才能创建的数据
        self.yt_st_num_list = yt_st_num_list
        self.yd_st_list = yd_st_list
        # 1, offset表示把m_ν'对应到m_ν的偏移量
        offset_dict_list = []
        for yd in self.yd_list:
            offset_dict = {}
            sum_num = 0
            bl_yd_list = self.bl_yd_list_dict.get(tuple(yd))
            for bl_yd in bl_yd_list:
                offset_dict[tuple(bl_yd)] = sum_num
                sum_num += yt_st_num_list[yd_st_list.index(bl_yd)]
            offset_dict_list.append(offset_dict)
        # 对于Sn=1，没有它。因为每次计算Sn和St的，所以S2时要保护一下
        self.offset_dict_list = offset_dict_list  # Sn=1, [{}]

    def get_yt_num(self, yd):
        return self.yt_num_list[self.yd_list.index(yd)]

    def get_yt_st_num(self, yd_st):
        return self.yt_st_num_list[self.yd_st_list.index(yd_st)]

    def get_eigenvalue(self, yd):
        return self.eigenvalue_list[self.yd_list.index(yd)]

    def get_tilde(self, yd):
        return self.yd_tilde_list[self.yd_list.index(yd)]

    def get_phase_factor_list(self, yd):
        return self.phase_factor_list_list[self.yd_list.index(yd)]

    def get_butler_phase_factor(self, yd, yd_st):
        return self.butler_phase_factor_dict_dict[tuple(yd)].get(tuple(yd_st), None)

    def get_offset(self, yd_st, yd):
        # 返回数字合法，返回None则不合法
        return self.offset_dict_list[self.yd_list.index(yd)].get(tuple(yd_st))

    def quick_calc_m(self, m_st, yd_st, yd):
        # 不是从yt开始数m，而是通过分支律，快速计算m
        # 是core.young_tableaux.py: quickly_calc_young_table_in_decreasing_page_order的特定情况简化
        return self.get_offset(yd_st, yd) + m_st

    def get_bl_yds(self, yd):
        return self.bl_yd_list_dict[tuple(yd)]

    # def get_bl_rows(self, yd):
    #     return self.bl_row_list_dict[tuple(yd)]

    def get_d3_k4_by_ϵ_key(self, ϵ_key: str):
        return self.ϵ_key2groups_dict[ϵ_key]

    def get_ϵ_key_by_d3_k4(self, d3: tuple, k4: tuple):
        return self.groups2ϵ_key_dict[(d3, k4)]

    def get_ν_combination_of_meta_σμν_tuple_list_by_σμ(self, σ, μ):
        if (σ, μ) in self.σμ_combination_of_meta_σμν_tuple_list:
            return self.ν_combination_of_meta_σμν_tuple_list_list[
                self.σμ_combination_of_meta_σμν_tuple_list.index((σ, μ))]
        else:
            return []


class CalcHelper(object):
    """ISFHelper, CGCHelper, EHelper公有的函数"""

    def __init__(self):
        self.s_n = None
        self.s_t = None
        # 按照ϵ的生成顺序，将σμν组合分成6组。第一个参数是σμν组合的字符串表示；第二个参数是σμν的顺序，第三个参数是该顺序下，每个变量是否取～
        self.σμν_0 = (("σμν", group_d3[0], group_k4[0]), ("σ~μ~ν", group_d3[0], group_k4[1]),
                      ("μσν", group_d3[1], group_k4[0]), ("μ~σ~ν", group_d3[1], group_k4[1]))  # m_ν=1
        self.σμν_1 = (("σ~μν~", group_d3[0], group_k4[2]), ("σμ~ν~", group_d3[0], group_k4[3]),
                      ("μ~σν~", group_d3[1], group_k4[2]), ("μσ~ν~", group_d3[1], group_k4[3]))  # m_ν=h_ν
        self.σμν_2 = (("σνμ", group_d3[3], group_k4[0]), ("σ~ν~μ", group_d3[3], group_k4[1]),
                      ("νσμ", group_d3[4], group_k4[0]), ("ν~σ~μ", group_d3[4], group_k4[1]))  # m_μ=1
        self.σμν_3 = (("σ~νμ~", group_d3[3], group_k4[2]), ("σν~μ~", group_d3[3], group_k4[3]),
                      ("ν~σμ~", group_d3[4], group_k4[2]), ("νσ~μ~", group_d3[4], group_k4[3]))  # m_μ=h_μ
        self.σμν_4 = (("νμσ", group_d3[2], group_k4[0]), ("ν~μ~σ", group_d3[2], group_k4[1]),
                      ("μνσ", group_d3[5], group_k4[0]), ("μ~ν~σ", group_d3[5], group_k4[1]))  # m_σ=1
        self.σμν_5 = (("ν~μσ~", group_d3[2], group_k4[2]), ("νμ~σ~", group_d3[2], group_k4[3]),
                      ("μ~νσ~", group_d3[5], group_k4[2]), ("μν~σ~", group_d3[5], group_k4[3]))  # m_σ=h_σ

    def enable_now_s_n(self, s_n):
        if not isinstance(s_n, int) or s_n < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_n={} must be int and >= {}".format(s_n, min(min_s_n_of_isf, min_s_n_of_cgc)))
        self.s_n = s_n
        self.s_t = s_n - 1


class ISFHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(ISFHelper, self).__init__()

    def calc_meta_isf_dict_include_save(self, data_σ_μ, ν_st, row_index_list, data_sn, data_st, old_isf_info_dict):
        """计算ISF中，属于元σμν的部分

        对于包含元σμν的组合σμ，计算其中元组合那部分的元ISF（无法通过对称性生成）
        使用如下公式：
        （主公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-192
        （相位公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-195

        正则ISF索引的全部参数为：σ σ' μ μ' ν τ ν' τ'
        表示：|σ σ'> * |μ μ'> 的结果中，|ντ ν'τ'>，的ISF系数平方

        返回flag, meta_isf_square_dict:
        meta_isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                                "cols":[[ν], ([ν], τ), ...],  # 有自由度tuple，无自由度list
                                "isf": isf_square_matrix}  # sp.Matrix([len(rows), len(cols)])，且len(rows)=len(cols)
        """
        meta_isf_start_time = time.time()
        # 只象征性检查Sn、St
        if not all(self.s_n == s_n for s_n in [data_sn.s_n, data_σ_μ.s_n]) \
                or not all(self.s_t == s_t for s_t in [data_st.s_n, data_σ_μ.s_t]):
            err_msg = "input wrong data with s_n={} should eq all(data_sn.s_n={}, data_σ_μ.s_n={}) " \
                      "and s_t={} should eq all(data_st.s_n={}, data_σ_μ.s_t={}) " \
                      "but not".format(self.s_n, data_sn.s_n, data_σ_μ.s_n, self.s_t, data_st.s_n, data_σ_μ.s_t)
            logger.error(err_msg)
            return False, err_msg

        if ν_st == [1]:  # 初始情况
            isf_2_square_dict = {"rows": [([1], [1])],
                                 "cols": [[2]] if data_σ_μ.σ == data_σ_μ.μ else [[1, 1]],
                                 "isf": sp.Matrix([[1]])}
            finish_cols = [[2]] if data_σ_μ.σ == data_σ_μ.μ else [[1, 1]]
            data_σ_μ.register_real_calc_isf_ν_st_list(ν_st)
            data_σ_μ.register_isf_square_dict_and_finish_cols_by_ν_st(ν_st, isf_2_square_dict, finish_cols)
            _, msg = save_isf(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st, isf_2_square_dict, 0)
            return True, None

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("@@@@ catch isf [5, 1, 1] [4, 2, 1] [3, 1, 1, 1]")

        start_time_4_193a = time.time()
        flag, isf_matrix_with_4_193a = self._calc_isf_matrix_4_193a(row_index_list, ν_st, data_sn, data_σ_μ, data_st)
        if not flag:
            err_msg = "_calc_isf_matrix_4_193a by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
                      "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix_with_4_193a)
            logger.error(err_msg)
            return False, err_msg
        speed_time_4_193a = int(time.time() - start_time_4_193a)
        isf_matrix = isf_matrix_with_4_193a
        # logger.debug("@@isf_matrix={} for data_σ_μ={}, ν_st={} is isf_matrix_with_4_193a"
        #              "".format(isf_matrix, data_σ_μ, ν_st))

        start_time_last_m = time.time()
        flag, isf_matrix_with_last_m = \
            self._calc_isf_matrix_with_last_m(row_index_list, ν_st, data_sn, data_σ_μ, data_st)
        if not flag:
            err_msg = "_calc_isf_matrix_with_last_m by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
                      "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix_with_last_m)
            logger.error(err_msg)
            return False, err_msg
        speed_time_last_m = int(time.time() - start_time_last_m)
        if isf_matrix != isf_matrix_with_last_m:
            logger.error("@@isf_matrix={} != isf_matrix_with_last_m={} with data_σ_μ={}, ν_st={}"
                         "".format(isf_matrix, isf_matrix_with_last_m, data_σ_μ, ν_st))

        start_time_first_m = time.time()
        flag, isf_matrix_with_first_m = self._calc_isf_matrix_with_first_m(row_index_list, ν_st, data_sn, data_σ_μ)
        if not flag:
            err_msg = "_calc_isf_matrix_with_first_m by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
                      "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix_with_first_m)
            logger.error(err_msg)
            return False, err_msg
        speed_time_first_m = int(time.time() - start_time_first_m)
        if isf_matrix != isf_matrix_with_first_m:
            logger.error("isf_matrix={} != isf_matrix_with_first_m={} with data_σ_μ={}, ν_st={}"
                         "".format(isf_matrix, isf_matrix_with_first_m, data_σ_μ, ν_st))

        if self.s_n >= 6 and speed_time_4_193a != 0 and speed_time_last_m != 0 and speed_time_first_m != 0:
            logger.debug("isf_matrix_with_4_193a using time={}\n"
                         "isf_matrix_with_last_m using time={}\n"
                         "isf_matrix_with_first_m using time={}\n"
                         "".format(speed_time_4_193a, speed_time_last_m, speed_time_first_m))

        # flag, isf_matrix_with_last_m = \
        #     self._calc_isf_matrix_with_last_m(row_index_list, ν_st, data_sn, data_σ_μ, data_st)
        # if not flag:
        #     err_msg = "_calc_isf_matrix_with_last_m by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
        #               "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix_with_last_m)
        #     logger.error(err_msg)
        #     return False, err_msg
        # if isf_matrix != isf_matrix_with_last_m:
        #     logger.error("isf_matrix={} != isf_matrix_with_last_m={} with data_σ_μ={}, ν_st={}"
        #                  "".format(isf_matrix, isf_matrix_with_last_m, data_σ_μ, ν_st))

        # flag, isf_matrix_with_4_193a = self._calc_isf_matrix_4_193a(row_index_list, ν_st, data_sn, data_σ_μ, data_st)
        # if not flag:
        #     err_msg = "_calc_isf_matrix_4_193a by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
        #               "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix)
        #     logger.error(err_msg)
        #     return False, err_msg
        # if isf_matrix != isf_matrix_with_4_193a:
        #     logger.error("isf_matrix={} != isf_matrix_with_4_193a={} with data_σ_μ={}, ν_st={}"
        #                  "".format(isf_matrix, isf_matrix_with_4_193a, data_σ_μ, ν_st))

        # flag, isf_matrix_with_4_193d = self._calc_isf_matrix_4_193d(row_index_list, ν_st, data_sn, data_σ_μ, data_st)
        # if not flag:
        #     err_msg = "isf_matrix_with_4_193d by row_index_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
        #               "with msg={}".format(row_index_list, ν_st, data_sn, data_σ_μ, isf_matrix)
        #     logger.error(err_msg)
        #     return False, err_msg
        # if isf_matrix != isf_matrix_with_4_193d:
        #     logger.error("isf_matrix={} != isf_matrix_with_4_193d={} with data_σ_μ={}, ν_st={}"
        #                  "".format(isf_matrix, isf_matrix_with_4_193d, data_σ_μ, ν_st))

        # 由于本征值是已知的，直接去算nullspace更简单。所以下面的代码注释掉
        # try:
        #     # Return list of triples (eigenval, multiplicity, eigenspace) eigenval升序
        #     '''
        #     p.s.
        #     Matrix([[1, 0, 2], [0, 3, 0], [2, 0, 1]]).eigenvects() return
        #     [(-1, 1, [Matrix([[-1], [ 0], [ 1]])]),
        #      (3,  2, [Matrix([[0],  [1],  [0]]),
        #               Matrix([[1],  [0],  [1]])])]
        #     注意这里的结果：单根是正交未归一的；多根不一定正交，也不归一
        #     '''
        #     # 到这里isf_matrix是可以通过符号计算，得到无误差的矩阵的。
        #     # 但后面，求本征值的时候，理论上不一定还能保持符号计算了（一元五次方程无公式解）
        #     # 不用怕，我们这里还有个终极大招，就是用np数值化计算理论上必然是整数的本征值，再回来使用符号计算
        #     # 或者，直接用cg_series预先拿到eigenvalues，只计算需要的eigenvectors就可以了!!!
        #     eigen_tuple_list = isf_matrix.eigenvects()
        # except Exception as e:
        #     logger.error(Exception(e))
        #     err_msg = "calc eigenvects meet fail with isf_matrix={} s_n={}".format(isf_matrix, self.s_n)
        #     logger.error(err_msg)
        #     return False, err_msg

        λ_ν_st = data_st.get_eigenvalue(ν_st)

        # 下面的代码改成先有col了，用λ(col)点名eigen_tuple_list，就可以不解矩阵，直接用mat+λ求eigenvectors
        if old_isf_info_dict is not False:
            meta_isf_square_dict = copy.deepcopy(old_isf_info_dict["data"])  # 这里浪费一点点，可以确保入参不会被改
            col_index_list = meta_isf_square_dict.get("cols")
            finish_cols = copy.deepcopy(old_isf_info_dict["flags"]["finish_cols"])
        else:
            col_index_list = self.calc_col_indexes(data_σ_μ, ν_st, data_sn)
            meta_isf_square_dict = {"rows": row_index_list,
                                    "cols": col_index_list,
                                    "isf": sp.zeros(len(row_index_list))}
            finish_cols = []
        ν_flag = None
        for tmp_col in col_index_list:
            if tmp_col in finish_cols:  # 算过了不用算（理论上应该不会出现）
                continue
            # tmp_col is [ν] or ([ν], τ)，但是，它不是最终的col，因为我们要重新指定τ
            ν = tmp_col if isinstance(tmp_col, list) else tmp_col[0]
            if ν == ν_flag:  # 相同的ν，不同的τ，只需做一次，就可以得到所有
                continue
            else:
                ν_flag = ν
            if ν not in data_σ_μ.meta_ν_list_of_σμ:  # 本函数只计算元σμν，非元的跳过不管（后面专门有只计算对称的）
                continue
                
            λ_ν = data_sn.get_eigenvalue(ν)
            e_value = λ_ν - λ_ν_st  # 这里就不用怕同λ的ν冲突了，因为col已经把不合理ν排除了
            e_vectors = self._calc_eigenspace_by_eigenvalue(isf_matrix, e_value)
            τ_max = len(e_vectors)
            # (e_value, τ_max, e_vectors) = list(filter(lambda x: x[0] == e_value, eigen_tuple_list))[0]
            if τ_max > 1:
                logger.debug("τ_max={} case with σ_μ_ν=({}{}{}), ν_st={} "
                             "\nand e_value={}, e_vectors={}, λ_ν={} - λ_ν_st={}"
                             "".format(τ_max, data_σ_μ.σ, data_σ_μ.μ, ν, ν_st, e_value, e_vectors, λ_ν, λ_ν_st))

            # 这里是用[::-1]倒序而不是正序就是为了对齐书中给的多重根自由度选择
            # TODO 可惜的是，就算如此，都无法保证soe_vectors与前人在相位上一致，原因有2：
            # 1，多自由度情况下，前人以什么顺序开展正交化是不确定的
            # 2，前人的位相可能是它的整体旋转，这样，任何一个e_vectors中的分量，都不在结果中
            soe_vectors = self._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if τ_max > 1 \
                else [self._calc_orthogonalization_vector(e_vectors[0])]
            if τ_max > 1:
                logger.debug("soe_vectors={}".format(soe_vectors))

            # 注意，τ_tmp_list, isf_phase_vector_list不仅可以带None占位，它还可以不是升序的，因为它对应着soe_vectors顺序

            # logger.debug("data_σ_μ={}, ν={}, ν_st={}, row_index_list={}, τ_max={}"
            #              "".format(data_σ_μ, ν, ν_st, row_index_list, τ_max))

            flag, τ_tmp_list, isf_phase_vector_list = \
                self._calc_isf_τ_and_phase_vector_list(soe_vectors, ν, ν_st, row_index_list, τ_max,
                                                       data_sn, data_st, data_σ_μ)
            if not flag:
                err_msg = "_calc_isf_τ_and_phase_vector_list meet error " \
                          "with soe_vectors={}, ν={}, ν_st={}, row_index_list={}, τ_max={} " \
                          "data_sn={}, data_st={}, data_σ_μ={}, " \
                          "msg={}".format(soe_vectors, ν, ν_st, row_index_list, τ_max,
                                          data_sn, data_st, data_σ_μ, τ_tmp_list)
                logger.error(err_msg)
                return False, err_msg
            if τ_max > 1:
                logger.debug("τ_tmp_list={}, isf_phase_vector_list={}".format(τ_tmp_list, isf_phase_vector_list))

            # 按照升序填空
            for τ_tmp in range(1, τ_max + 1):  # 因为按照soe_vectors对应的τ乱序了，所以需要正序
                τ_tmp = None if τ_max == 1 else τ_tmp
                τ_tmp_index = τ_tmp_list.index(τ_tmp)
                phase_vector = isf_phase_vector_list[τ_tmp_index]
                # 计算单列的ISF
                isf_square = sp.Matrix([sp.sign(i) * i**2 for i in phase_vector])
                real_col = ν if τ_max == 1 else (ν, τ_tmp)  # 区别于tmp_col，它是真的了
                real_col_index = col_index_list.index(real_col)
                meta_isf_square_dict["isf"][:, real_col_index] = isf_square
                # 填finish_cols
                finish_cols.append(real_col)  # 它的顺序在这里不重要

        # 注册元ISF信息和数据备用
        data_σ_μ.register_real_calc_isf_ν_st_list(ν_st)
        data_σ_μ.register_isf_square_dict_and_finish_cols_by_ν_st(ν_st, meta_isf_square_dict, finish_cols)

        # 保存 # 注意：这里的ISF仍然可能是不完整的
        meta_isf_speed_time = int(time.time() - meta_isf_start_time)
        old_adding_time = 0 if old_isf_info_dict is False else old_isf_info_dict["flags"]["speed_time"]
        if len(finish_cols) == len(row_index_list):  # 当完整的时候，它们长度也相等
            flag, msg = save_isf(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st, meta_isf_square_dict,
                                 meta_isf_speed_time + old_adding_time)
            if not flag:
                err_msg = "save meta isf meet error with s_n={}, meta_σ={}, meta_μ={}, ν_st={}, " \
                          "meta_isf_square_dict={}, msg={}" \
                          "".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st, meta_isf_square_dict, msg)
                logger.error(err_msg)
                return False, err_msg
        else:
            flag, msg = save_isf(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st, meta_isf_square_dict,
                                 meta_isf_speed_time + old_adding_time, finish_cols=finish_cols)
            if not flag:
                err_msg = "save meta isf meet error with s_n={}, meta_σ={}, meta_μ={}, ν_st={}, " \
                          "sym_isf_square_dict={}, finish_cols={}, msg={}" \
                          "".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st, meta_isf_square_dict, finish_cols, msg)
                logger.error(err_msg)
                return False, err_msg

        return True, None

    @staticmethod
    def _calc_eigenspace_by_eigenvalue(mat, eigenvalue):
        """
        从sympy源码摘抄并简化的一段，用来计算eigenspace
        Get a basis for the eigenspace for a particular eigenvalue
        """
        m = mat - sp.eye(mat.rows) * eigenvalue
        ret = m.nullspace()
        # the nullspace for a real eigenvalue should be non-trivial.
        if len(ret) == 0:
            raise NotImplementedError(
                "Can't evaluate eigenvector for eigenvalue %s" % eigenvalue)
        return ret

    def _calc_isf_τ_and_phase_vector_list(self, soe_vectors, ν, ν_st, row_index_list, τ_max,
                                          data_sn, data_st, data_σ_μ):
        """ 将经过施密特正交归一化手续后的本征矢量调整为Yamanouchi相位，同时确定非第一分支的τ
        返回结果是按照soe_vectors对应排序的

        对于abs对应不上的，需要调整到与第一分支对齐的soe_vector

        1，将分支律的第一分支首个非零系数调整为正（绝对相位）
        2，非第一分支的，参考第一分支做相对调整（相对相位）
        """
        τ_tmp_list = []  # 这个tmp是相对于完全完成的τ是没有None的，但它int部分的τ就是最终确定的τ
        phase_vector_list = []

        if data_sn.get_bl_yds(ν)[0] == ν_st:  # ν_st击中ν的第一分支

            # if _debug_condition(data_σ_μ, ν_st):
            #     logger.warning("@@@@ 绝对相位 ν={} 的第一分支是 ν_st={}".format(ν, ν_st))

            # 当ν_st为ν第一分支的时候，令首个非零系数符号为+，称为绝对相位
            for τ_tmp, soe_vector in zip(range(1, τ_max + 1), soe_vectors):
                τ_tmp = None if τ_max == 1 else τ_tmp
                _, no_0_element = self._get_first_no_0_number_from_vector(soe_vector)
                if no_0_element is None:
                    err_msg = "find schmidt_orthogonalization_vector={} all approximately equal to 0, pls check".format(
                        soe_vector)
                    logger.error(err_msg)
                    return False, err_msg, None
                τ_tmp_list.append(τ_tmp)
                absolute_phase = sp.sign(no_0_element)
                phase_vector_list.append(absolute_phase * soe_vector)

            return True, τ_tmp_list, phase_vector_list

        else:  # 未击中情况，则相对相位需要参考它的分支律第一分支（ν_st_fbl）的绝对相位
            # 被参考的第一分支被称为fbl: first branching law
            '''
            这里的计算思路是，根据式4-195c，既可以用第一分支，计算其他分支的isf；也可以用其他分支，反向计算第一分支
            所以，如果仅仅差一个相位的isf已经得到，就可以用它反推第一分支，相同则不改相位；不同则*-1
            若不仅是符号不同，则使用正推，调整到符合首分支的相位
            另外，非第一分支的τ也可以用此方法确定（也是反推，直到确定属于哪个τ）
            '''
            ν_st_fbl = data_sn.get_bl_yds(ν)[0]  # ν_st_fbl就是ν_st对应ν的第一分支，它们σ、μ虽然不同，但σ'μ'ν相同

            # if _debug_condition(data_σ_μ, ν_st):
            #     logger.warning("@@@@ 相对相位 ν={} 的第一分支不是 ν_st={} 而是 ν_st_fbl={}".format(ν, ν_st, ν_st_fbl))

            # 三个m的意义分别是：
            # TODO 优化 或者解释一下它的本质（利用ν_st_st本质上是去掉Sn和Sn-1，重新生长）
            # 1 m_ν：               令当前ν_st（非ν的第一分支）的m'取1时，ν的m；
            # 2 m_ν_st_fbl：        令当前ν_st的第一分支ν_st_st的m''取1时，ν_st_fbl（ν的第一分支）的m'
            #                       (ν_st_st是ν_st的第一分支，可以证明它也是ν_st_fbl的分支，但不一定是第一分支)
            # 3,m_ν_by_m_ν_st_fbl： 按照2中方法，确定m_ν_st_fbl后（也意味着确定了ν_st_fbl的杨盘），
            #                       ν对应的m
            m_ν_st = 1
            m_ν = data_sn.quick_calc_m(m_ν_st, ν_st, ν)
            ν_st_st = data_st.get_bl_yds(ν_st)[0]
            m_ν_st_fbl = data_st.quick_calc_m(m_ν_st, ν_st_st, ν_st_fbl)
            m_ν_by_m_ν_st_fbl = data_sn.quick_calc_m(m_ν_st_fbl, ν_st_fbl, ν)

            # 因为按照ν_st的Yamanouchi序开启循环，所以前面的fbl分支，必然先于ν_st被计算
            flag, isf_fbl_dict = load_isf(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st_fbl,
                                          is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_isf meet error with self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, ν_st_fbl={}," \
                          "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st_fbl, isf_fbl_dict)
                logger.error(err_msg)
                return False, err_msg, None

            isf_fbl_row_list, isf_fbl_col_list, isf_fbl_square = \
                isf_fbl_dict["rows"], isf_fbl_dict["cols"], isf_fbl_dict["isf"]

            # if _debug_condition(data_σ_μ, ν_st):
            #     logger.warning("@@@@ isf_fbl_row_list={}, isf_fbl_col_list={}, isf_fbl_square={} with "
            #                    "data_σ_μ.σ={}, data_σ_μ.μ={}, ν_st_fbl={}".format(
            #         isf_fbl_row_list, isf_fbl_col_list, isf_fbl_square, data_σ_μ.σ, data_σ_μ.μ, ν_st_fbl))

            # 将第一分支(fbl)的τ分配给当前分支，并充分利用结果，把能确定的phase定下来
            '''
            为了区分τ，我们采取一个最笨的办法，就是：
            1，维护一个flag_uncertain_τ_num，用来判断是否完成目标
            2，维护τ_by_soe_vectors，用来表示当前每个vector的τ范围（循环结束后应该是len1的）
            3，维护phase_by_soe_vectors，用来充分利用结果，把已经计算出的phase写进去，后面就不用再算了
            for 循环所有soe_vectors的行（fbl作为答案，行号和soe_vectors的是一致的）
                3.1，判断flag_uncertain_τ_num
                3.2，把当前未确定τ的列 对应的isf_fbl_square提取出来
                3.3，计算abs集合，判断是否需要跳过该结果
                
                for 循环soe_vectors的所有列（未完全确定τ的列，它们有的有可能已经被缩小范围了）
                    3.3.0，跳过已经确定τ的列
                    3.3.3.1，用整体符号待定的isf计算对应行的isf_fbl
                    3.3.3.2，缩小范围
                    3.3.3.3，判断是否全部τ都已经被区分
                        3-1，区分，退出
                        3-2，未区分，开始下一轮循环
            4，对未确定phase的列，补充phase
            '''
            τ_fbl_list = [None] if τ_max == 1 else list(range(1, τ_max + 1))  # load的isf必然按照τ升序排列
            # 初始所有τ_fbl_list都是候选，随着逐渐缩小范围直至全部确定  # 它与soe_vectors顺序是对应的
            τ_candidate = [set(τ_fbl_list) for i in range(τ_max)]  # 它用列表解析就是深拷贝，用*就是浅拷贝
            phase_candidate = [None] * τ_max  # 确定一个填一个，允许循环结束后有未确定的
            fbl_cols_useful = [ν] if τ_max == 1 else list((ν, i) for i in τ_fbl_list)  # 指的是fbl中能用于判断τ和phase的那部分
            fbl_cols_useful_index = [isf_fbl_col_list.index(col) for col in fbl_cols_useful]
            τ_make_sure_set = {None} if τ_max == 1 else set()  # 也作为统一旗标，避免指标太多，逻辑混乱
            flag_soe = True  # 标志本soe是否与首分支相位一致，一致，则使用；否则，要用首分支计算全部phase_vector

            for row_index, row in enumerate(isf_fbl_row_list):  # 既是fbl的，也是soe_vectors的
                if len(τ_make_sure_set) == τ_max:
                    break
                isf_fbl_square_by_row = isf_fbl_square[row_index, :]
                isf_fbl_square_useful_list = [isf_fbl_square_by_row[i] for i in fbl_cols_useful_index]
                isf_fbl_square_abs_τ_set_dict = {}  # abs(isf_fbl_square): set(τ_possible)
                for single_fbl_τ, single_isf_fbl_square in zip(τ_fbl_list, isf_fbl_square_useful_list):
                    if single_fbl_τ in τ_make_sure_set:
                        continue
                    isf_fbl_square_abs = abs(single_isf_fbl_square)
                    if isf_fbl_square_abs in isf_fbl_square_abs_τ_set_dict:
                        isf_fbl_square_abs_τ_set_dict[isf_fbl_square_abs] &= {single_fbl_τ}
                    else:
                        isf_fbl_square_abs_τ_set_dict[isf_fbl_square_abs] = {single_fbl_τ}
                isf_fbl_square_abs_diff_set = set(isf_fbl_square_abs_τ_set_dict.keys())
                if len(isf_fbl_square_abs_diff_set) == 1:
                    # 绝对值全相等，说明无区分τ的能力，跳过
                    # TODO 如果全部跳过，说明τ任意，程序暂时未考虑出现此情况
                    # TODO 还有一种未考虑情况如：a1=[1/2, -1/2] a2=[-1/2, 1/2], a3=[1/3, -2/3], a4=[2/3, -1/3]...
                    continue

                for single_τ_candidate_set, (col_index, soe_vector) in zip(τ_candidate, enumerate(soe_vectors)):
                    if len(single_τ_candidate_set) == 1:  # 跳过确定好的τ
                        continue
                    # 用当前待定τ的列soe_vector，去计算已知的isf_fbl_square
                    flag, isf_fbl_another = \
                        self._calc_isf_by_known_isf(soe_vector, row_index_list, ν_st, m_ν_st, m_ν,
                                                    row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl, ν, data_sn, data_σ_μ)
                    if not flag:
                        err_msg = "calc isf_fbl_another meet error " \
                                  "with soe_vector={}, row_index_list={}, ν_st={}, m_ν_st={}, m_ν={}, " \
                                  "row={}, ν_st_fbl={}, m_ν_st_fbl={}, m_ν_by_m_ν_st_fbl={}, " \
                                  "ν={}, data_sn={}, data_σ_μ={}, " \
                                  "msg={}".format(soe_vector, row_index_list, ν_st, m_ν_st, m_ν,
                                                  row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                  ν, data_sn, data_σ_μ, isf_fbl_another)
                        logger.error(err_msg)
                        return False, err_msg, None
                    isf_fbl_another_square_abs = isf_fbl_another ** 2
                    # abs检查
                    if isf_fbl_another_square_abs not in isf_fbl_square_abs_diff_set:
                        if τ_max == 1:  # 无自由度的情况下，要求必须相等，所以此时要报错
                            err_msg = "isf_fbl_another_square_abs={} not in isf_fbl_square_abs_diff_set={}, pls check" \
                                      "with isf_fbl_another={}, row={}, soe_vector={}, data_sn={}, data_st={}, " \
                                      "data_σ_μ={}".format(isf_fbl_another_square_abs, isf_fbl_square_abs_diff_set,
                                                           isf_fbl_another, row, soe_vector, data_sn, data_st, data_σ_μ)
                            logger.error(err_msg)
                            return False, err_msg, None
                        else:  # 有自由度情况，遇到一次不相等，则推翻这次的自由度，采用4-195c正算
                            # msg = "isf_fbl_another_square_abs={} not in isf_fbl_square_abs_diff_set={}", \
                            #       "with σ_μ_ν=({}{}{}), ν_st={}, and isf_fbl_another={}, row={} " \
                            #       "using soe_vector={}" \
                            #       "".format(isf_fbl_another_square_abs, isf_fbl_square_abs_diff_set,
                            #                 data_σ_μ.σ, data_σ_μ.μ, ν, ν_st, isf_fbl_another, row, soe_vector)
                            # logger.debug(msg)
                            logger.debug("will break and use 4-195c to calc(1)")
                            flag_soe = False
                            break
                    # 缩小τ范围
                    isf_fbl_square_abs_τ_set = isf_fbl_square_abs_τ_set_dict[isf_fbl_another_square_abs]
                    single_τ_candidate_set &= isf_fbl_square_abs_τ_set  # 用旧范围交集新范围
                    if len(single_τ_candidate_set) == 1:
                        # 完全确定了一个τ，flag_uncertain_τ_num和isf_fbl_square_abs_τ_set_dict都要随之改变
                        isf_fbl_square_abs_τ_set_dict[isf_fbl_another_square_abs] -= single_τ_candidate_set
                        τ_make_sure_set |= single_τ_candidate_set
                    # 顺带确定phase（只有完全确定了τ的才能定phase）
                    if isf_fbl_another_square_abs != 0 and len(single_τ_candidate_set) == 1:
                        # phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(soe_vector[row_index])
                        fbl_col = ν if τ_max == 1 else (ν, list(single_τ_candidate_set)[0])
                        fbl_col_index = isf_fbl_col_list.index(fbl_col)
                        phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(isf_fbl_square[row_index, fbl_col_index])
                        if phase_candidate[col_index] is None:
                            phase_candidate[col_index] = phase_tmp_new
                        else:
                            if phase_candidate[col_index] != phase_tmp_new:  # 既然算了，不检查一下浪费了
                                err_msg = "find phase_tmp_new={} not eq old={} " \
                                          "with data_sn={}, data_st={}, data_σ_μ={}".format(
                                    phase_tmp_new, phase_candidate[col_index], data_sn, data_st, data_σ_μ)
                                logger.error(err_msg)
                                return False, err_msg, None

            # if _debug_condition(data_σ_μ, ν_st):
            #     logger.warning("@@@@ τ_candidate={}".format(τ_candidate))
            #     logger.warning("@@@@ phase_candidate={}".format(phase_candidate))
            #     logger.warning("@@@@ soe_vectors={}".format(soe_vectors))

            if flag_soe:
                # 到这里，soe_vectors的τ已经完全确定，phase部分确定，下面补充未确定的phase就可以了
                for τ_tmp_set, phase_tmp, soe_vector in zip(τ_candidate, phase_candidate, soe_vectors):
                    # 检查τ
                    if len(τ_tmp_set) != 1:
                        err_msg = "find τ_tmp_set={} not union with data_sn={}, data_st={}, data_σ_μ={}".format(
                            τ_tmp_set, data_sn, data_st, data_σ_μ)
                        logger.error(err_msg)
                        return False, err_msg, None

                        '''
                         core.isf_and_cgc:isf_and_cgc.py:259 calc_meta_isf_dict_include_save meet error by data_σ_μ=σ=[5, 2], μ=[4, 2, 1] data class, ν_st=[3, 1, 1, 1], isf_row_index_list=[([5, 1], [4, 1, 1]), ([5, 1], [3, 2, 1]), ([4, 2], [4, 2]), ([4, 2], [4, 1, 1]), ([4, 2], [3, 2, 1], 1), ([4, 2], [3, 2, 1], 2)], data_si=s_n=7 data class, data_st=s_n=6 data class, old_isf_info_dict={'data': {'rows': [([5, 1], [4, 1, 1]), ([5, 1], [3, 2, 1]), ([4, 2], [4, 2]), ([4, 2], [4, 1, 1]), ([4, 2], [3, 2, 1], 1), ([4, 2], [3, 2, 1], 2)], 'cols': [([4, 1, 1, 1], 1), ([4, 1, 1, 1], 2), ([3, 2, 1, 1], 1), ([3, 2, 1, 1], 2), ([3, 2, 1, 1], 3), [3, 1, 1, 1, 1]], 'isf': Matrix([
                        [0, 0, 0, 0, 0,  5/14],
                        [0, 0, 0, 0, 0, -5/28],
                        [0, 0, 0, 0, 0, 27/70],
                        [0, 0, 0, 0, 0,  1/14],
                        [0, 0, 0, 0, 0,     0],
                        [0, 0, 0, 0, 0, 1/140]])}, 'flags': {'speed_time': 0, 'finish_cols': [[3, 1, 1, 1, 1]]}} with msg=_calc_isf_τ_and_phase_vector_list meet error with soe_vectors=[Matrix([
                        [            0],
                        [         -1/2],
                        [   -sqrt(6)/5],
                        [2*sqrt(10)/15],
                        [   sqrt(10)/6],
                        [        -7/30]]), Matrix([
                        [          0],
                        [  sqrt(6)/4],
                        [       1/10],
                        [ sqrt(15)/5],
                        [          0],
                        [-sqrt(6)/20]]), Matrix([
                        [ sqrt(3)/3],
                        [sqrt(6)/12],
                        [      -1/2],
                        [         0],
                        [         0],
                        [ sqrt(6)/4]])], ν=[3, 2, 1, 1], ν_st=[3, 1, 1, 1], row_index_list=[([5, 1], [4, 1, 1]), ([5, 1], [3, 2, 1]), ([4, 2], [4, 2]), ([4, 2], [4, 1, 1]), ([4, 2], [3, 2, 1], 1), ([4, 2], [3, 2, 1], 2)], τ_max=3 data_sn=s_n=7 data class, data_st=s_n=6 data class, data_σ_μ=σ=[5, 2], μ=[4, 2, 1] data class, msg=find τ_tmp_set=set() not union with data_sn=s_n=7 data class, data_st=s_n=6 data class, data_σ_μ=σ=[5, 2], μ=[4, 2, 1] data class
                        '''

                    τ_tmp = list(τ_tmp_set)[0]
                    τ_tmp_list.append(τ_tmp)
                    # 补漏phase
                    if phase_tmp is not None:
                        phase_vector_list.append(phase_tmp * soe_vector)
                    else:
                        fbl_col = ν if τ_max == 1 else (ν, τ_tmp)
                        fbl_col_index = isf_fbl_col_list.index(fbl_col)
                        isf_fbl_square_same_τ = isf_fbl_square[:, fbl_col_index]
                        first_no_0_isf_fbl_row_index, first_no_0_isf_fbl_square = \
                            self._get_first_no_0_number_from_vector(isf_fbl_square_same_τ)

                        # if _debug_condition(data_σ_μ, ν_st):
                        #     logger.warning("@@@@ first_no_0_isf_fbl_row_index={}".format(first_no_0_isf_fbl_row_index))
                        #     logger.warning("@@@@ first_no_0_isf_fbl_square={}".format(first_no_0_isf_fbl_square))

                        first_no_0_isf_fbl_row = isf_fbl_row_list[first_no_0_isf_fbl_row_index]
                        flag, isf_fbl_another = \
                            self._calc_isf_by_known_isf(soe_vector, row_index_list, ν_st, m_ν_st, m_ν,
                                                        first_no_0_isf_fbl_row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                        ν, data_sn, data_σ_μ)
                        if not flag:
                            err_msg = "complete isf_fbl_another meet error " \
                                      "with soe_vector={}, row_index_list={}, ν_st={}, m_ν_st={}, m_ν={}, " \
                                      "first_no_0_isf_fbl_row={}, ν_st_fbl={}, m_ν_st_fbl={}, m_ν_by_m_ν_st_fbl={}, " \
                                      "ν={}, data_sn={}, data_σ_μ={}, " \
                                      "msg={}".format(soe_vector, row_index_list, ν_st, m_ν_st, m_ν,
                                                      first_no_0_isf_fbl_row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                      ν, data_sn, data_σ_μ, isf_fbl_another)
                            logger.error(err_msg)
                            return False, err_msg, None
                        isf_fbl_another_square_abs = isf_fbl_another ** 2
                        # abs检查
                        if isf_fbl_another_square_abs != abs(first_no_0_isf_fbl_square):
                            if τ_max == 1:  # 无自由度的情况下，要求必须相等，所以此时要报错
                                err_msg = "isf_fbl_another_square_abs={} not eq abs(first_no_0_isf_fbl_square)={} " \
                                          "with isf_fbl_another={}, first_no_0_isf_fbl_row={}, soe_vector={}, " \
                                          "data_sn={}, data_st={}, data_σ_μ={}".format(
                                    isf_fbl_another_square_abs, first_no_0_isf_fbl_square,
                                    isf_fbl_another, first_no_0_isf_fbl_row, soe_vector, data_sn, data_st, data_σ_μ)
                                logger.error(err_msg)
                                return False, err_msg, None
                            else:  # 有自由度情况，遇到一次不相等，则推翻这次的自由度，采用4-195c正算
                                msg = "isf_fbl_another_square_abs={} not eq abs(first_no_0_isf_fbl_square)={}", \
                                      "with σ_μ_ν=({}{}{}), ν_st={} and isf_fbl_another={}, first_no_0_isf_fbl_row={}" \
                                      "using soe_vector={}" \
                                      "".format(isf_fbl_another_square_abs, first_no_0_isf_fbl_square,
                                                data_σ_μ.σ, data_σ_μ.μ, ν, ν_st,
                                                isf_fbl_another, first_no_0_isf_fbl_row, soe_vector)
                                logger.debug(msg)
                                logger.debug("will break and use 4-195c to calc(2)")
                                flag_soe = False
                                break
                        # TODO 这里应该是拿待定相位的τ算出来的isf_fbl_another和确定的isf_fbl_another比较啊
                        # phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(soe_vector[first_no_0_isf_fbl_row_index])
                        phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(first_no_0_isf_fbl_square)
                        phase_vector_list.append(phase_tmp_new * soe_vector)

            if not flag_soe:
                τ_tmp_list = []
                phase_vector_list = []

                # logger.warning("@@@@ renew isf for ν={}, ν_st={}, data_sn={}, data_st={}, data_σ_μ={}".format(
                #     ν, ν_st, data_sn, data_st, data_σ_μ))

                # 使用4-195c重新正算phase_vector，放弃seo_vector
                for τ_tmp in range(1, τ_max + 1):  # 此时，列按照fbl的就是正序了，且它必须有非零τ
                    τ_tmp_list.append(τ_tmp)
                    col = (ν, τ_tmp)
                    isf_fbl_square_by_col = isf_fbl_square[:, isf_fbl_col_list.index(col)]
                    isf_fbl_col_vector = sp.Matrix([sp.sign(i) * sp.sqrt(abs(i)) for i in isf_fbl_square_by_col])
                    phase_vector_tmp = []
                    for row in row_index_list:
                        flag, single_phase_vector_new = \
                            self._calc_isf_by_known_isf(isf_fbl_col_vector, isf_fbl_row_list, ν_st_fbl,
                                                        m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                        row, ν_st, m_ν_st, m_ν,
                                                        ν, data_sn, data_σ_μ)
                        if not flag:
                            err_msg = "calc single_phase_vector_new meet error " \
                                      "with isf_fbl_col_vector={}, isf_fbl_row_list={}, ν_st_fbl={}, " \
                                      "m_ν_st_fbl={}, m_ν_by_m_ν_st_fbl={} " \
                                      "row={}, ν_st={}, m_ν_st={}, m_ν={}, " \
                                      "ν={}, data_sn={}, data_σ_μ={}, " \
                                      "msg={}".format(isf_fbl_col_vector, isf_fbl_row_list, ν_st_fbl,
                                                      m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                      row, ν_st, m_ν_st, m_ν,
                                                      ν, data_sn, data_σ_μ, single_phase_vector_new)
                            logger.error(err_msg)
                            return False, err_msg, None
                        phase_vector_tmp.append(single_phase_vector_new)
                    phase_vector_list.append(sp.Matrix(phase_vector_tmp))

                # logger.warning("@@@@ renew isf phase_vector_list={}".format(phase_vector_list))
                # logger.warning("@@@@ old soe_vectors={}".format(soe_vectors))

                # TODO delete it(因为这部分回推，对计算是多余的，它只是初次写这段代码的检查)
                # 利用new的，再回推fbl，反正以下自洽
                for τ_tmp in range(1, τ_max + 1):
                    col = (ν, τ_tmp)
                    first_isf_fbl_row = isf_fbl_row_list[0]
                    first_isf_fbl_square = isf_fbl_square[:, isf_fbl_col_list.index(col)][0]
                    phase_vector = phase_vector_list[τ_tmp - 1]
                    flag, isf_fbl_another = \
                        self._calc_isf_by_known_isf(phase_vector, row_index_list, ν_st, m_ν_st, m_ν,
                                                    first_isf_fbl_row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                    ν, data_sn, data_σ_μ)
                    if sp.sign(isf_fbl_another) * isf_fbl_another**2 != first_isf_fbl_square:
                        logger.error("$$$$ 回推不成立！")
                        err_msg = "first_isf_fbl_square={}, isf_fbl_another={}".format(
                            first_isf_fbl_square, isf_fbl_another)
                        logger.error(err_msg)
                logger.warning("$$$$ 回推成立！with σ_μ_ν=({}{}{}), ν_st={}".format(data_σ_μ.σ, data_σ_μ.μ, ν, ν_st))

        return True, τ_tmp_list, phase_vector_list

    def _calc_isf_by_known_isf(self, k_isf_vector, k_isf_rows, k_ν_st, k_m_ν_st, k_m_ν,
                               u_isf_row, u_ν_st, u_m_ν_st, u_m_ν_by_k_m_ν_st,
                               ν, data_sn, data_σ_μ):
        """
        根据式4-195c(英文书4-195b)，根据已知ISF计算未知

        其中，正算为：
        由确定相位的第一分支(k_isf)全部，逐个直接推导准确的其余分支(u_isf)
        反算为：
        使用待定的非第一分支(k_isf)全部，逐个回推确定的第一分支(u_isf)，
        如果abs相等，说明只需要调整整体符号；如果abs不等，则说明需要调整符号和整体旋转

        输入：
        第一行：known独有参数
        第二行：u独有参数
        第三行：公有参数
        """
        # TODO 应该也可以用4-193a的思想优化
        u_σ_st, u_μ_st, u_τ_st = u_isf_row if len(u_isf_row) == 3 else (u_isf_row[0], u_isf_row[1], None)
        flag, u_cgc_square_st_dict = load_cgc(self.s_t, u_σ_st, u_μ_st, u_ν_st, u_τ_st, u_m_ν_st,
                                              is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get u_cgc_square_st_dict with s_t={}, u_σ_st={}, u_μ_st={}, u_ν_st={}, u_τ_st={}, " \
                      "u_m_ν_st={} meet error with " \
                      "msg={}".format(self.s_t, u_σ_st, u_μ_st, u_ν_st, u_τ_st, u_m_ν_st, u_cgc_square_st_dict)
            logger.error(err_msg)
            return False, err_msg

        sum_3_loop = 0
        offset_σ_u = data_sn.get_offset(u_σ_st, data_σ_μ.σ)
        offset_μ_u = data_sn.get_offset(u_μ_st, data_σ_μ.μ)
        for k_isf_row, k_vector_element in zip(k_isf_rows, k_isf_vector):
            k_σ_st, k_μ_st, k_τ_st = k_isf_row if len(k_isf_row) == 3 else (k_isf_row[0], k_isf_row[1], None)
            flag, k_cgc_square_st_dict = load_cgc(self.s_t, k_σ_st, k_μ_st, k_ν_st, k_τ_st, k_m_ν_st,
                                                  is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get k_cgc_square_st_dict with s_t={}, k_σ_st={}, k_μ_st={}, k_ν_st={}, k_τ_st={}, " \
                          "k_m_ν_st={} meet error " \
                          "with msg={}".format(self.s_t, k_σ_st, k_μ_st, k_ν_st, k_τ_st, k_m_ν_st, k_cgc_square_st_dict)
                logger.error(err_msg)
                return False, err_msg
            sum_3_loop_part = 0  # 暂时先不优化此处，因为cgc的储存也可能被优化
            offset_σ_k = data_sn.get_offset(k_σ_st, data_σ_μ.σ)
            offset_μ_k = data_sn.get_offset(k_μ_st, data_σ_μ.μ)
            for (u_m_σ_st, u_m_μ_st), u_cgc_square_st_element in u_cgc_square_st_dict.items():
                u_m_σ = offset_σ_u + u_m_σ_st
                u_m_μ = offset_μ_u + u_m_μ_st
                for (k_m_σ_st, k_m_μ_st), k_cgc_square_st_element in k_cgc_square_st_dict.items():
                    k_m_σ = offset_σ_k + k_m_σ_st
                    k_m_μ = offset_μ_k + k_m_μ_st
                    in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][u_m_σ - 1, k_m_σ - 1]
                    in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][u_m_μ - 1, k_m_μ - 1]
                    sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
                                       * sp.sign(k_cgc_square_st_element) * sp.sign(u_cgc_square_st_element) \
                                       * sp.sqrt(abs(k_cgc_square_st_element * u_cgc_square_st_element))
            sum_3_loop += sum_3_loop_part * k_vector_element

        flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in")  # (Sn-1, Sn)的对换
        if not flag:
            err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
                      "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
            logger.error(err_msg)
            return False, err_msg
        # in_matrix_ν_m_m_element = in_matrix_ν[k_m_ν - 1, u_m_ν_by_k_m_ν_st - 1]
        in_matrix_ν_m_m_element = in_matrix_ν[u_m_ν_by_k_m_ν_st - 1, k_m_ν - 1]  # 这个改动其实没变化，因为in矩阵是对称的
        u_isf_element = sum_3_loop / in_matrix_ν_m_m_element
        return True, u_isf_element

    @staticmethod
    def _get_first_no_0_number_from_vector(vector):
        """提取vector中首个非0的数字和index"""
        for i, number in enumerate(vector):
            if number != 0:
                return i, number
        return None, None

    @staticmethod
    def _calc_orthogonalization_vector(vector):
        sum_square = sum(i**2 for i in vector)
        return vector / sp.sqrt(sum_square)

    @staticmethod
    def _calc_schmidt_orthogonalization_tricky(multi_vectors):
        """使用GramSchmidt得到正交归一的多重根本征矢量
        这里，为了使存在自由度的结果对齐《群表示论的新途径》，选用一个tricky的技巧:
        method1：就是取分母之和表示复杂度，返回复杂度最小的一组
        method2：取分母平方和
        method3：取分子分母各自的平方和
        """
        # soe_1 = sp.GramSchmidt(multi_vectors, True)
        # soe_2 = sp.GramSchmidt(multi_vectors[::-1], True)
        # soe_1_denominator, soe_2_denominator = 0, 0
        # for mat_1, mat_2 in zip(soe_1, soe_2):
        #     soe_1_denominator += sum(i.as_numer_denom()[1] for i in mat_1)  # 0的'分母'会被取为1，不影响计算啦
        #     soe_2_denominator += sum(i.as_numer_denom()[1] for i in mat_2)
        # rst = soe_1 if soe_1_denominator <= soe_2_denominator else soe_2
        #
        # if rst == soe_1:
        #     logger.warning("kill soe_2={} with soe_1_denominator={} <= soe_2_denominator={}"
        #                    "".format(soe_2, soe_1_denominator, soe_2_denominator))
        # else:
        #     logger.warning("kill soe_1={} with soe_1_denominator={} > soe_2_denominator={}"
        #                    "".format(soe_1, soe_1_denominator, soe_2_denominator))
        # return rst

        # soe_default = sp.GramSchmidt(multi_vectors, True)
        # soe_default_denominator = sum(sum(i.as_numer_denom()[1] for i in single_vector)
        #                               for single_vector in soe_default)  # 0的'分母'都会被取为1，所以不影响比较
        # if len(multi_vectors) > 2:
        #     logger.warning("%%%%%%%% " + "{}".format(len(multi_vectors)) * 10)
        soe = None
        soe_denominator = None
        soe_permutation = None
        diff_kill = []
        same_kill = []

        vectors_index_permutations = list(permutations(range(len(multi_vectors)), len(multi_vectors)))
        for vector_index in vectors_index_permutations:
            vector_candidate = [multi_vectors[i] for i in vector_index]
            soe_candidate = sp.GramSchmidt(vector_candidate, True)
            soe_candidate_denominator = sum(sum(i.as_numer_denom()[1] for i in single_vector)
                                            for single_vector in soe_candidate)  # 0的'分母'都会被取为1，所以不影响比较

            if soe_denominator is None:  # 初始化
                soe = soe_candidate
                soe_denominator = soe_candidate_denominator
                soe_permutation = vector_index

            elif soe_candidate_denominator < soe_denominator:  # 替换  # 还要把当前的same_kill剪切到diff_kill中
                # logger.warning("replace larger_soe={} with old_denominator={} "
                #                "to smaller_soe={} with new_denominator={} by method1, "
                #                "the new vector_index={} "
                #                "".format(soe, soe_denominator, soe_candidate, soe_candidate_denominator,
                #                          soe_permutation))
                # 先剪切
                diff_kill.append(soe_denominator)
                for k in same_kill:
                    diff_kill.append(k[1])
                same_kill = []
                # 再把当前的放到kill_list里
                diff_kill.append(soe_candidate_denominator)  # 抛弃，并记录到diff_kill中
                # logger.warning("kill vector_candidate={} with {}>{} by method1"
                #                "".format(soe_candidate, soe_candidate_denominator, soe_denominator))
                # 最后替换
                soe_denominator = soe_candidate_denominator
                soe = soe_candidate
                soe_permutation = vector_index

            elif soe_candidate_denominator == soe_denominator:  # 抛弃，但记录到same_kill中
                # logger.warning("kill vector_candidate={} with {}=={} by method1"
                #                "".format(soe_candidate, soe_candidate_denominator, soe_denominator))
                same_kill.append([soe_candidate, soe_candidate_denominator, vector_index])
            else:
                diff_kill.append(soe_candidate_denominator)  # 抛弃，并记录到diff_kill中
                # logger.warning("kill vector_candidate={} with {}>{} by method1"
                #                "".format(soe_candidate, soe_candidate_denominator, soe_denominator))

            # 下面这段到else为止，是错误的，故意反过来判断
            # elif soe_candidate_denominator > soe_denominator:  # 替换  # 还要把当前的same_kill剪切到diff_kill中
            #     logger.warning("replace larger_soe={} with old_denominator={} "
            #                    "to smaller_soe={} with new_denominator={} by method1, "
            #                    "the new vector_index={} "
            #                    "".format(soe, soe_denominator, soe_candidate, soe_candidate_denominator,
            #                              soe_permutation))
            #     # 先剪切
            #     diff_kill.append(soe_denominator)
            #     for k in same_kill:
            #         diff_kill.append(k[1])
            #     same_kill = []
            #     # 后替换
            #     soe_denominator = soe_candidate_denominator
            #     soe = soe_candidate
            #     soe_permutation = vector_index
            #
            # elif soe_candidate_denominator == soe_denominator:  # 抛弃，但记录到same_kill中
            #     logger.warning("kill vector_candidate={} with {}=={} by method1"
            #                    "".format(soe_candidate, soe_candidate_denominator, soe_denominator))
            #     same_kill.append([soe_candidate, soe_candidate_denominator, vector_index])
            # else:
            #     diff_kill.append(soe_candidate_denominator)  # 抛弃，并记录到diff_kill中
            #     logger.warning("kill vector_candidate={} with {}<{} by method1"
            #                    "".format(soe_candidate, soe_candidate_denominator, soe_denominator))

        logger.warning("len of vectors_index_permutations={}".format(len(vectors_index_permutations)))
        logger.warning("diff kill len={}, min={} and larger_list={}"
                       "".format(len(diff_kill), soe_denominator, sorted(diff_kill)))
        logger.warning("same kill len={}, list={}".format(len(same_kill), same_kill))
        logger.warning("$$min soe={} with soe_denominator={}, soe_permutation={} by method1"
                       "".format(soe, soe_denominator, soe_permutation))

        return soe

    # @staticmethod
    # def _calc_ν_by_λ_and_bl(λ_ν, eigenvalue_list, yd_list, bl_yd_list_dict, ν_st):
    #     """根据本征值航道，得到ν"""
    #     count_λ_ν = eigenvalue_list.count(λ_ν)
    #     if count_λ_ν == 0:
    #         err_msg = "λ_ν={} should in eigenvalue_list={} but not, pls check".format(λ_ν, eigenvalue_list)
    #         logger.error(err_msg)
    #         raise Exception(err_msg)
    #     elif count_λ_ν == 1:
    #         return yd_list[eigenvalue_list.index(λ_ν)]
    #
    #     else:  # 多于1个，需要用"航道"确定
    #         rst_list = []
    #         for ev, yd_candidate in zip(eigenvalue_list, yd_list):
    #             if ev != λ_ν:
    #                 continue
    #             bl_yds = bl_yd_list_dict[tuple(yd_candidate)]
    #             if ν_st in bl_yds:  # 核心就是判断yd_candidate的分支律中是否含有ν_st
    #                 rst_list.append(yd_candidate)
    #         if len(rst_list) == 1:
    #             return rst_list[0]
    #         else:
    #             # TODO 这里先报错，暂不确定是否会出现这种情况
    #             msg = "get rst_list={} len > 1, pls check, " \
    #                   "with λ_ν={}, eigenvalue_list={}, yd_list={}, bl_yd_list_dict={}, ν_st={}".format(
    #                 rst_list, λ_ν, eigenvalue_list, yd_list, bl_yd_list_dict, ν_st)
    #             logger.warning(msg)
    #             raise Exception(msg)

    # @staticmethod
    # def calc_row_indexes_tmp(data_σ_μ, ν_st, data_st, is_tmp=True):
    #     """计算ISF表格的行的意义，它是的bl_σ, bl_μ, τ'的列表
    #     形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
    #     row_index_tmp_list = []  # [(bl_σ, bl_μ, τ'), ...]
    #     for bl_σ, bl_μ in product(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ):
    #         single_ν_cg_series = data_σ_μ.get_cg_series_st(bl_σ, bl_μ, ν_st, data_st.yd_list)
    #         if single_ν_cg_series == 0:
    #             continue
    #         part_rst_list = [(bl_σ, bl_μ, τ_st) if single_ν_cg_series >= 2 else
    #                          ((bl_σ, bl_μ, None) if is_tmp is True else (bl_σ, bl_μ))
    #                          for τ_st in range(1, single_ν_cg_series + 1)]  # [1, 2, ..., single_ν_cg_series]
    #         row_index_tmp_list += part_rst_list
    #     return row_index_tmp_list  # 当σ * μ中不包含ν_st时，结果为[]

    @staticmethod
    def calc_row_indexes(data_σ_μ, ν_st, data_st):
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, τ'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        row_index_list = []  # [(bl_σ, bl_μ, τ'), ...]
        for bl_σ, bl_μ in product(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ):
            single_ν_cg_series = data_σ_μ.get_cg_series_st(bl_σ, bl_μ, ν_st, data_st.yd_list)
            if single_ν_cg_series == 0:
                continue
            part_rst_list = [(bl_σ, bl_μ, τ_st) if single_ν_cg_series >= 2 else
                             (bl_σ, bl_μ) for τ_st in range(1, single_ν_cg_series + 1)]
            row_index_list += part_rst_list
        return row_index_list  # 当σ * μ中不包含ν_st时，结果为[]

    @staticmethod
    def calc_col_indexes(data_σ_μ, ν_st, data_sn):
        """计算ISF表格的列的意义，它是的ν, τ的列表
        形如[[4,1], ([3,2],1), ([3,2],2), [3,1,1], ...]"""
        col_index_list = []  # [(ν, τ), ...]
        for yd, cg_series in zip(data_sn.yd_list, data_σ_μ.cg_series_list):
            if cg_series == 0:
                continue
            yd_bl = data_sn.get_bl_yds(yd)
            if ν_st in yd_bl:
                if cg_series == 1:
                    col_index_list.append(yd)
                else:
                    for i in range(1, 1+cg_series):
                        col_index_list.append((yd, i))
        return col_index_list

    @lru_cache()
    def _load_cgc_with_m1_by_input_json(self, s_n, j_σ, j_μ, j_ν, τ):
        """lru版本的load_cgc，通过将入参更改为json.dumps实现序列化/json.loads还原"""
        σ, μ, ν = map(json.loads, [j_σ, j_μ, j_ν])
        return load_cgc_with_m1(s_n, σ, μ, ν, τ)

    def _load_cgc_with_m1_lru(self, s_n, σ, μ, ν, τ):
        j_σ, j_μ, j_ν = map(json.dumps, [σ, μ, ν])
        return self._load_cgc_with_m1_by_input_json(s_n, j_σ, j_μ, j_ν, τ)

    def _cgc_st_2_cgc_m_dict(self, cgc_st_tuple, data_sn, data_σ_μ):
        """取得cg_st_dict，并且，根据分支律，将m_st转化为对应的m
        注意：返回的字典，仅仅是把m_st按照分支律准化成了它的m对应，value仍然是st的cgc"""
        σ_st, μ_st, ν_st, τ_st, _ = cgc_st_tuple
        rst_dict = {}
        flag, cgc_st_square_dict = self._load_cgc_with_m1_lru(self.s_t, σ_st, μ_st, ν_st, τ_st)
        if not flag:
            err_msg = "load_cgc fail by s_t={}, σ_st={}, μ_st={}, ν_st={}, τ_st={} with msg={}".format(
                self.s_t, σ_st, μ_st, ν_st, τ_st, cgc_st_square_dict)
            logger.error(err_msg)
            return False, err_msg
        offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
        offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
        for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
            m_σ = offset_σ + m_σ_st
            m_μ = offset_μ + m_μ_st
            rst_dict[(m_σ, m_μ,)] = cgc_st_square

        return True, rst_dict

    # def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, ν_st=None):
    #     """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
    #     matrix_element = 0
    #     # 将St的cgc根据分支律上升到Sn，并整形成py序
    #     # 左
    #     _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_left, data_sn, data_σ_μ)
    #     left_tmp_dict = {}
    #     for left_key, factor_cgc_left in factor_cgc_left_dict.items():
    #         left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
    #         left_tmp_dict[left_key] = left_tmp
    #     # 右
    #     if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
    #         right_tmp_dict = left_tmp_dict
    #     else:
    #         _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_right, data_sn, data_σ_μ)
    #         right_tmp_dict = {}
    #         for right_key, factor_cgc_right in factor_cgc_right_dict.items():
    #             right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
    #             right_tmp_dict[right_key] = right_tmp
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("@@@@ left_tmp_dict={} with cgc_st_tuple_left={}, and factor_cgc_left_dict={}, "
    #                        "right_tmp_dict={} with cgc_st_tuple_right={}"
    #                        "".format(left_tmp_dict, cgc_st_tuple_left, factor_cgc_left_dict,
    #                                  right_tmp_dict, cgc_st_tuple_right))
    #
    #     # 计算matrix_element
    #     for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
    #         for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
    #             in_element_sum = 0
    #             for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
    #                 σ_in_element = data_σ_μ.in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
    #                 μ_in_element = data_σ_μ.in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1, m_μ_right - 1]
    #                 in_element_sum += σ_in_element * μ_in_element
    #
    #                 if _debug_condition(data_σ_μ, ν_st):
    #                     logger.warning("@@@@ ({},{}), m_σ_left={}, m_σ_right={}, σ_in_element={}, "
    #                                    "m_μ_left={}, m_μ_right={}, μ_in_element={}"
    #                                    "".format(i, self.s_n, m_σ_left, m_σ_right, σ_in_element,
    #                                              m_μ_left, m_μ_right, μ_in_element))
    #
    #             if _debug_condition(data_σ_μ, ν_st):
    #                 logger.warning("@@@@ in_element_sum={}, left_tmp={}, right_tmp={}"
    #                                "".format(in_element_sum, left_tmp, right_tmp))
    #
    #             matrix_element += in_element_sum * left_tmp * right_tmp
    #
    #     return True, matrix_element
    #
    # def _calc_isf_matrix(self, row_index_list, ν_st, data_sn, data_σ_μ):
    #     """计算ISF的本征矩阵
    #     """
    #     matrix_div = len(row_index_list)
    #     isf_matrix = sp.zeros(matrix_div)
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("_calc_isf_matrix")
    #         for i in range(1, self.s_n):
    #             logger.warning("σ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_σ_dict[(i, self.s_n)]))
    #             logger.warning("μ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_μ_dict[(i, self.s_n)]))
    #
    #     # 构建主对角元
    #     for i, row in enumerate(row_index_list):
    #         cgc_st_tuple = (row[0], row[1], ν_st, None, 1) if len(row) == 2 else (row[0], row[1], ν_st, row[2], 1)
    #         flag, isf_matrix_element = \
    #             self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, ν_st)
    #         if not flag:
    #             err_msg = "get isf_matrix_element by cgc_st_tuple={}, cgc_st_tuple={}, data_sn={}, data_σ_μ={}" \
    #                       "with msg={}".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
    #             logger.error(err_msg)
    #             return False, err_msg
    #
    #         if _debug_condition(data_σ_μ, ν_st):
    #             logger.warning("@@@@ i={}, i={}, isf_matrix_element={}"
    #                            "".format(i, i, isf_matrix_element))
    #
    #         isf_matrix[i, i] = isf_matrix_element
    #
    #     # 构建其他矩阵元素
    #     if matrix_div >= 2:
    #         for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
    #             cgc_st_tuple_left = (left_st[0], left_st[1], ν_st, None, 1) if len(left_st) == 2 \
    #                 else (left_st[0], left_st[1], ν_st, left_st[2], 1)
    #             cgc_st_tuple_right = (right_st[0], right_st[1], ν_st, None, 1) if len(right_st) == 2 \
    #                 else (right_st[0], right_st[1], ν_st, right_st[2], 1)
    #             flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple_left, cgc_st_tuple_right,
    #                                                                      data_sn, data_σ_μ, ν_st)
    #             if not flag:
    #                 err_msg = "get isf_matrix_element by cgc_st_tuple_left={}, cgc_st_tuple_right={}, " \
    #                           "data_sn={}, data_σ_μ={} " \
    #                           "with msg={}".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn,
    #                                                data_σ_μ, isf_matrix_element)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #             if _debug_condition(data_σ_μ, ν_st):
    #                 logger.warning("@@@@ left_i={}, right_i={}, isf_matrix_element={}"
    #                                "".format(left_i, right_i, isf_matrix_element))
    #
    #             isf_matrix[left_i, right_i] = isf_matrix_element
    #             isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭
    #
    #     return True, isf_matrix

    def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, ν_st=None, i_n=None):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_left, data_sn, data_σ_μ)
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
            left_tmp_dict[left_key] = left_tmp
        # 右
        if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
            right_tmp_dict = left_tmp_dict
        else:
            _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_right, data_sn, data_σ_μ)
            right_tmp_dict = {}
            for right_key, factor_cgc_right in factor_cgc_right_dict.items():
                right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
                right_tmp_dict[right_key] = right_tmp

        # if _debug_condition(data_σ_μ, ν_st):
        #     logger.warning("@@@@ left_tmp_dict={} with cgc_st_tuple_left={}, and factor_cgc_left_dict={}, "
        #                    "right_tmp_dict={} with cgc_st_tuple_right={}"
        #                    "".format(left_tmp_dict, cgc_st_tuple_left, factor_cgc_left_dict,
        #                              right_tmp_dict, cgc_st_tuple_right))

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                σ_in_element = data_σ_μ.in_matrix_σ_dict[i_n][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
                μ_in_element = data_σ_μ.in_matrix_μ_dict[i_n][m_μ_left - 1, m_μ_right - 1]
                matrix_element += σ_in_element * μ_in_element * left_tmp * right_tmp

        return True, matrix_element

    def _calc_isf_matrix_with_first_m(self, row_index_list, ν_st, data_sn, data_σ_μ):
        """计算ISF的本征矩阵
        """
        matrix_div = len(row_index_list)
        isf_matrix = sp.zeros(matrix_div)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("####_calc_isf_matrix_with_first_m####")
            # for i in range(1, self.s_n):
            #     logger.warning("σ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_σ_dict[(i, self.s_n)]))
            #     logger.warning("μ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_μ_dict[(i, self.s_n)]))

        for j in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
            i_n = (j, self.s_n)
            sub_isf_matrix = sp.zeros(matrix_div)

            # 构建主对角元
            for i, row in enumerate(row_index_list):
                cgc_st_tuple = (row[0], row[1], ν_st, None, 1) if len(row) == 2 else (row[0], row[1], ν_st, row[2], 1)
                flag, isf_matrix_element = \
                    self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, ν_st, i_n)
                if not flag:
                    err_msg = "get isf_matrix_element by cgc_st_tuple={}, cgc_st_tuple={}, data_sn={}, data_σ_μ={}" \
                              "with msg={}".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                sub_isf_matrix[i, i] = isf_matrix_element

            # 构建其他矩阵元素
            if matrix_div >= 2:
                for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
                    cgc_st_tuple_left = (left_st[0], left_st[1], ν_st, None, 1) if len(left_st) == 2 \
                        else (left_st[0], left_st[1], ν_st, left_st[2], 1)
                    cgc_st_tuple_right = (right_st[0], right_st[1], ν_st, None, 1) if len(right_st) == 2 \
                        else (right_st[0], right_st[1], ν_st, right_st[2], 1)
                    flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple_left, cgc_st_tuple_right,
                                                                             data_sn, data_σ_μ, ν_st, i_n)
                    if not flag:
                        err_msg = "get isf_matrix_element by cgc_st_tuple_left={}, cgc_st_tuple_right={}, " \
                                  "data_sn={}, data_σ_μ={} " \
                                  "with msg={}".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn,
                                                       data_σ_μ, isf_matrix_element)
                        logger.error(err_msg)
                        return False, err_msg
                    sub_isf_matrix[left_i, right_i] = isf_matrix_element
                    sub_isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ i_n={}, sub_isf_matrix={}"
                               "".format(i_n, sub_isf_matrix))

            isf_matrix += sub_isf_matrix

        return True, isf_matrix

    def _cgc_st_2_cgc_m_dict_with_point_m(self, cgc_st_tuple, data_sn, data_σ_μ):
        """取得cg_st_dict，并且，根据分支律，将m_st转化为对应的m
        注意：返回的字典，仅仅是把m_st按照分支律准化成了它的m对应，value仍然是st的cgc"""
        σ_st, μ_st, ν_st, τ_st, m_st = cgc_st_tuple
        rst_dict = {}
        flag, cgc_st_square_dict = load_cgc(self.s_t, *cgc_st_tuple, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "load_cgc fail by s_t={}, cgc_st_tuple={} with msg={}".format(
                self.s_t, cgc_st_tuple, cgc_st_square_dict)
            logger.error(err_msg)
            return False, err_msg
        offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
        offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
        for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
            m_σ = offset_σ + m_σ_st
            m_μ = offset_μ + m_μ_st
            rst_dict[(m_σ, m_μ,)] = cgc_st_square

        return True, rst_dict

    # def _calc_isf_matrix_element_with_last_m(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, ν_st=None):
    #     """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
    #     matrix_element = 0
    #     # 将St的cgc根据分支律上升到Sn，并整形成py序
    #     # 左
    #     _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_left, data_sn, data_σ_μ)
    #     left_tmp_dict = {}
    #     for left_key, factor_cgc_left in factor_cgc_left_dict.items():
    #         left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
    #         left_tmp_dict[left_key] = left_tmp
    #     # 右
    #     if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
    #         right_tmp_dict = left_tmp_dict
    #     else:
    #         _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_right, data_sn, data_σ_μ)
    #         right_tmp_dict = {}
    #         for right_key, factor_cgc_right in factor_cgc_right_dict.items():
    #             right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
    #             right_tmp_dict[right_key] = right_tmp
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("@@@@ left_tmp_dict={} with cgc_st_tuple_left={}, and factor_cgc_left_dict={}, "
    #                        "right_tmp_dict={} with cgc_st_tuple_right={}"
    #                        "".format(left_tmp_dict, cgc_st_tuple_left, factor_cgc_left_dict,
    #                                  right_tmp_dict, cgc_st_tuple_right))
    #
    #     # 计算matrix_element
    #     for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
    #         for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
    #             in_element_sum = 0
    #             for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
    #                 σ_in_element = data_σ_μ.in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
    #                 μ_in_element = data_σ_μ.in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1, m_μ_right - 1]
    #                 in_element_sum += σ_in_element * μ_in_element
    #
    #                 if _debug_condition(data_σ_μ, ν_st):
    #                     logger.warning("@@@@ ({},{}), m_σ_left={}, m_σ_right={}, σ_in_element={}, "
    #                                    "m_μ_left={}, m_μ_right={}, μ_in_element={}"
    #                                    "".format(i, self.s_n, m_σ_left, m_σ_right, σ_in_element,
    #                                              m_μ_left, m_μ_right, μ_in_element))
    #
    #             if _debug_condition(data_σ_μ, ν_st):
    #                 logger.warning("@@@@ in_element_sum={}, left_tmp={}, right_tmp={}"
    #                                "".format(in_element_sum, left_tmp, right_tmp))
    #
    #             matrix_element += in_element_sum * left_tmp * right_tmp
    #
    #     return True, matrix_element

    # def _calc_isf_matrix_with_last_m(self, row_index_list, ν_st, data_sn, data_σ_μ, data_st):
    #     """计算ISF的本征矩阵
    #     """
    #     matrix_div = len(row_index_list)
    #     isf_matrix = sp.zeros(matrix_div)
    #     h_ν_st = data_st.get_yt_num(ν_st)
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("#### _calc_isf_matrix_with_last_m")
    #         for i in range(1, self.s_n):
    #             logger.warning("σ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_σ_dict[(i, self.s_n)]))
    #             logger.warning("μ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_μ_dict[(i, self.s_n)]))
    #
    #     # 构建主对角元
    #     for i, row in enumerate(row_index_list):
    #         cgc_st_tuple = (row[0], row[1], ν_st, None, h_ν_st) if len(row) == 2 else \
    #             (row[0], row[1], ν_st, row[2], h_ν_st)
    #         flag, isf_matrix_element = self._calc_isf_matrix_element_with_last_m(cgc_st_tuple, cgc_st_tuple,
    #                                                                              data_sn, data_σ_μ, ν_st)
    #         if not flag:
    #             err_msg = "get _calc_isf_matrix_element_with_last_m by cgc_st_tuple={}, cgc_st_tuple={}, " \
    #                       "data_sn={}, data_σ_μ={} with msg={}" \
    #                       "".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
    #             logger.error(err_msg)
    #             return False, err_msg
    #
    #         if _debug_condition(data_σ_μ, ν_st):
    #             logger.warning("@@@@ i={}, i={}, isf_matrix_element={}".format(i, i, isf_matrix_element))
    #
    #         isf_matrix[i, i] = isf_matrix_element
    #
    #     # 构建其他矩阵元素
    #     if matrix_div >= 2:
    #         for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
    #             cgc_st_tuple_left = (left_st[0], left_st[1], ν_st, None, h_ν_st) if len(left_st) == 2 \
    #                 else (left_st[0], left_st[1], ν_st, left_st[2], h_ν_st)
    #             cgc_st_tuple_right = (right_st[0], right_st[1], ν_st, None, h_ν_st) if len(right_st) == 2 \
    #                 else (right_st[0], right_st[1], ν_st, right_st[2], h_ν_st)
    #             flag, isf_matrix_element = \
    #                 self._calc_isf_matrix_element_with_last_m(cgc_st_tuple_left, cgc_st_tuple_right,
    #                                                           data_sn, data_σ_μ, ν_st)
    #             if not flag:
    #                 err_msg = "get _calc_isf_matrix_element_with_last_m by cgc_st_tuple_left={}, " \
    #                           "cgc_st_tuple_right={}, data_sn={}, data_σ_μ={} with msg={}" \
    #                           "".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, isf_matrix_element)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #             if _debug_condition(data_σ_μ, ν_st):
    #                 logger.warning("@@@@ left_i={}, right_i={}, isf_matrix_element={}"
    #                                "".format(left_i, right_i, isf_matrix_element))
    #
    #             isf_matrix[left_i, right_i] = isf_matrix_element
    #             isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭
    #
    #     return True, isf_matrix
    #
    def _calc_isf_matrix_element_with_last_m(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ,
                                             ν_st=None, i_n=None):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_left, data_sn, data_σ_μ)
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
            left_tmp_dict[left_key] = left_tmp
        # 右
        if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
            right_tmp_dict = left_tmp_dict
        else:
            _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_right, data_sn, data_σ_μ)
            right_tmp_dict = {}
            for right_key, factor_cgc_right in factor_cgc_right_dict.items():
                right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
                right_tmp_dict[right_key] = right_tmp

        # if _debug_condition(data_σ_μ, ν_st):
        #     logger.warning("@@@@ left_tmp_dict={} with cgc_st_tuple_left={}, and factor_cgc_left_dict={}, "
        #                    "right_tmp_dict={} with cgc_st_tuple_right={}"
        #                    "".format(left_tmp_dict, cgc_st_tuple_left, factor_cgc_left_dict,
        #                              right_tmp_dict, cgc_st_tuple_right))

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                σ_in_element = data_σ_μ.in_matrix_σ_dict[i_n][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
                μ_in_element = data_σ_μ.in_matrix_μ_dict[i_n][m_μ_left - 1, m_μ_right - 1]
                matrix_element += σ_in_element * μ_in_element * left_tmp * right_tmp

        return True, matrix_element

    def _calc_isf_matrix_with_last_m(self, row_index_list, ν_st, data_sn, data_σ_μ, data_st):
        """计算ISF的本征矩阵
        """
        matrix_div = len(row_index_list)
        isf_matrix = sp.zeros(matrix_div)
        h_ν_st = data_st.get_yt_num(ν_st)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("####_calc_isf_matrix_with_last_m####")
            # for i in range(1, self.s_n):
            #     logger.warning("σ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_σ_dict[(i, self.s_n)]))
            #     logger.warning("μ_({},{})={}".format(i, self.s_n, data_σ_μ.in_matrix_μ_dict[(i, self.s_n)]))

        for j in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
            i_n = (j, self.s_n)
            sub_isf_matrix = sp.zeros(matrix_div)

            # 构建主对角元
            for i, row in enumerate(row_index_list):
                cgc_st_tuple = (row[0], row[1], ν_st, None, h_ν_st) if len(row) == 2 else \
                    (row[0], row[1], ν_st, row[2], h_ν_st)
                flag, isf_matrix_element = self._calc_isf_matrix_element_with_last_m(cgc_st_tuple, cgc_st_tuple,
                                                                                     data_sn, data_σ_μ, ν_st, i_n)
                if not flag:
                    err_msg = "get _calc_isf_matrix_element_with_last_m by cgc_st_tuple={}, cgc_st_tuple={}, " \
                              "data_sn={}, data_σ_μ={} with msg={}" \
                              "".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                sub_isf_matrix[i, i] = isf_matrix_element

            # 构建其他矩阵元素
            if matrix_div >= 2:
                for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
                    cgc_st_tuple_left = (left_st[0], left_st[1], ν_st, None, h_ν_st) if len(left_st) == 2 \
                        else (left_st[0], left_st[1], ν_st, left_st[2], h_ν_st)
                    cgc_st_tuple_right = (right_st[0], right_st[1], ν_st, None, h_ν_st) if len(right_st) == 2 \
                        else (right_st[0], right_st[1], ν_st, right_st[2], h_ν_st)
                    flag, isf_matrix_element = \
                        self._calc_isf_matrix_element_with_last_m(cgc_st_tuple_left, cgc_st_tuple_right,
                                                                  data_sn, data_σ_μ, ν_st, i_n)
                    if not flag:
                        err_msg = "get _calc_isf_matrix_element_with_last_m by cgc_st_tuple_left={}, " \
                                  "cgc_st_tuple_right={}, data_sn={}, data_σ_μ={} with msg={}" \
                                  "".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, isf_matrix_element)
                        logger.error(err_msg)
                        return False, err_msg
                    sub_isf_matrix[left_i, right_i] = isf_matrix_element
                    sub_isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ i_n={}, sub_isf_matrix={}"
                               "".format(i_n, sub_isf_matrix))

            isf_matrix += sub_isf_matrix

        return True, isf_matrix

    def _calc_sub_isf_matrix_element_4_193a(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ, ν_st=None):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化 选用4-193a优化
        4-193a的思想是每种m''得到一种值，所以每个m''派一个代表就可以了。它不需要对in求和了"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_left, data_sn, data_σ_μ)
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
            left_tmp_dict[left_key] = left_tmp
        # 右
        if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
            right_tmp_dict = left_tmp_dict
        else:
            _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict_with_point_m(cgc_st_tuple_right, data_sn, data_σ_μ)
            right_tmp_dict = {}
            for right_key, factor_cgc_right in factor_cgc_right_dict.items():
                right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
                right_tmp_dict[right_key] = right_tmp

        # if _debug_condition(data_σ_μ, ν_st):
        #     logger.warning("@@@@ left_tmp_dict={} with cgc_st_tuple_left={}, and factor_cgc_left_dict={}, "
        #                    "right_tmp_dict={} with cgc_st_tuple_right={}"
        #                    "".format(left_tmp_dict, cgc_st_tuple_left, factor_cgc_left_dict,
        #                              right_tmp_dict, cgc_st_tuple_right))

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                σ_in_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_left - 1, m_σ_right - 1]
                μ_in_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_left - 1, m_μ_right - 1]
                matrix_element += σ_in_element * μ_in_element * left_tmp * right_tmp

                # if _debug_condition(data_σ_μ, ν_st):
                #     logger.warning("@@@@ σ_in_element={}, μ_in_element={}, left_tmp={}, right_tmp={}"
                #                    "".format(σ_in_element, μ_in_element, left_tmp, right_tmp))

        return True, matrix_element

    def _calc_isf_matrix_4_193a(self, row_index_list, ν_st, data_sn, data_σ_μ, data_st):
        """计算ISF的本征矩阵 选用4-193a优化
        4-193a的思想是每种m''得到一种值，所以每个m''派一个代表就可以了"""
        if self.s_n <= 2:
            return self._calc_isf_matrix_with_first_m(row_index_list, ν_st, data_sn, data_σ_μ)

        matrix_div = len(row_index_list)
        isf_matrix = sp.zeros(matrix_div)
        h_ν_st = data_st.get_yt_num(ν_st)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("####_calc_isf_matrix_4_193a####")
            # logger.warning("σ_({},{})={}".format(self.s_t, self.s_n, data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)]))
            # logger.warning("μ_({},{})={}".format(self.s_t, self.s_n, data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)]))

        # ν_st的分支律
        ν_st_st_list = data_st.get_bl_yds(ν_st)
        for ν_st_st in ν_st_st_list:
            m_ν_st_st = 1
            m_ν_st = data_st.get_offset(ν_st_st, ν_st) + m_ν_st_st
            h_ν_st_st = data_st.get_yt_st_num(ν_st_st)
            sub_isf_matrix = sp.zeros(matrix_div)

            # 构建主对角元
            for i, row in enumerate(row_index_list):
                cgc_st_tuple = (row[0], row[1], ν_st, None, m_ν_st) if len(row) == 2 else \
                    (row[0], row[1], ν_st, row[2], m_ν_st)
                flag, isf_matrix_element = \
                    self._calc_sub_isf_matrix_element_4_193a(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, ν_st)
                if not flag:
                    err_msg = "_calc_sub_isf_matrix_element_4_193a by cgc_st_tuple={}, cgc_st_tuple={}, " \
                              "data_sn={}, data_σ_μ={} with msg={}" \
                              "".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                sub_isf_matrix[i, i] = isf_matrix_element

            # 构建其他矩阵元素
            if matrix_div >= 2:
                for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
                    cgc_st_tuple_left = (left_st[0], left_st[1], ν_st, None, m_ν_st) if len(left_st) == 2 \
                        else (left_st[0], left_st[1], ν_st, left_st[2], m_ν_st)
                    cgc_st_tuple_right = (right_st[0], right_st[1], ν_st, None, m_ν_st) if len(right_st) == 2 \
                        else (right_st[0], right_st[1], ν_st, right_st[2], m_ν_st)
                    flag, isf_matrix_element = \
                        self._calc_sub_isf_matrix_element_4_193a(cgc_st_tuple_left, cgc_st_tuple_right,
                                                                 data_sn, data_σ_μ, ν_st)
                    if not flag:
                        err_msg = "_calc_sub_isf_matrix_element_4_193a by cgc_st_tuple_left={}, " \
                                  "cgc_st_tuple_right={}, data_sn={}, data_σ_μ={} with msg={}" \
                                  "".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ,
                                            isf_matrix_element)
                        logger.error(err_msg)
                        return False, err_msg
                    sub_isf_matrix[left_i, right_i] = isf_matrix_element
                    sub_isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭

                    if _debug_condition(data_σ_μ, ν_st) and left_i == 3 and right_i == 5:
                        logger.debug("cgc_st_tuple_left={}, cgc_st_tuple_right={}"
                                     "".format(cgc_st_tuple_left, cgc_st_tuple_right))
                        flag, isf_matrix_element_anti = \
                            self._calc_sub_isf_matrix_element_4_193a(cgc_st_tuple_right, cgc_st_tuple_left,
                                                                     data_sn, data_σ_μ, ν_st)
                        logger.debug("isf_matrix_element = {} and isf_matrix_element_anti = {}"
                                     "".format(isf_matrix_element, isf_matrix_element_anti))
                        pass

            # 融合
            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ ν_st_st={}, h_ν_st_st={}, sub_isf_matrix={}"
                               "".format(ν_st_st, h_ν_st_st, sub_isf_matrix))

            isf_matrix += h_ν_st_st * sub_isf_matrix

        isf_matrix = isf_matrix * (self.s_n - 1) / h_ν_st

        return True, isf_matrix

    @staticmethod
    def _get_m_range_by_σμ_st_and_σμ_st_st(σμ_st_st, σμ_st_left, σμ_st_right, σμ, data_sn, data_st):
        """拿m''的范围去造m的范围（注意，是喂给range的范围），一次出左右两个m（左右st_st相同，st不同）"""
        last_m_st_st = data_st.get_yt_st_num(σμ_st_st)
        left_m_offset = data_sn.get_offset(σμ_st_left, σμ) + data_st.get_offset(σμ_st_st, σμ_st_left)
        right_m_offset = data_sn.get_offset(σμ_st_right, σμ) + data_st.get_offset(σμ_st_st, σμ_st_right)
        return (left_m_offset + 1, left_m_offset + last_m_st_st + 1), \
               (right_m_offset + 1, right_m_offset + last_m_st_st + 1)

    def _calc_sub_isf_matrix_element_4_193d(self, isf_st_tuple_left, isf_st_tuple_right, data_sn, data_st, data_σ_μ):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化 选用4-193d优化
        4-193d的思想是用ISF代替CGC去构建矩阵"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        _, factor_isf_st_left_dict = load_isf(self.s_t, *isf_st_tuple_left, is_flag_true_if_not_s_n=False)
        isf_st_left_cols = factor_isf_st_left_dict["cols"]
        isf_st_left_rows = factor_isf_st_left_dict["rows"]
        isf_st_left_square = factor_isf_st_left_dict["isf"]
        # 右
        if isf_st_tuple_left == isf_st_tuple_right:  # 如果相等，则直接使用left结果，无需重复计算
            factor_isf_st_right_dict = factor_isf_st_left_dict
        else:
            _, factor_isf_st_right_dict = load_isf(self.s_t, *isf_st_tuple_right, is_flag_true_if_not_s_n=False)
        isf_st_right_cols = factor_isf_st_right_dict["cols"]
        isf_st_right_rows = factor_isf_st_right_dict["rows"]
        isf_st_right_square = factor_isf_st_right_dict["isf"]

        # 计算matrix_element
        for left_i, left_row in enumerate(isf_st_left_rows):
            σ_st_st_left, μ_st_st_left = left_row[0], left_row[1]  # 暂时不需要τ_st_st_left
            for right_i, right_row in enumerate(isf_st_right_rows):
                σ_st_st_right, μ_st_st_right = right_row[0], right_row[1]  # 暂时不需要τ_st_st_left
                # 必须相等的σ''、μ''、ν''、τ''
                if σ_st_st_left != σ_st_st_right or μ_st_st_left != μ_st_st_right:
                    continue  # 用循环看起来会比较对称，其实点名也可以

                # logger.warning("#### σ_st_st_left={}, isf_st_tuple_left[0]={}, isf_st_tuple_right[0]={}, "
                #                "data_σ_μ.σ={}, data_sn={}, data_st={}"
                #                "".format(σ_st_st_left, isf_st_tuple_left[0], isf_st_tuple_right[0],
                #                          data_σ_μ.σ, data_sn, data_st))
                if _debug_condition(data_σ_μ):
                    logger.warning("σ_st_st_left,μ_st_st_left={}{}, | σ_st_st_right,μ_st_st_right={}{}"
                                   "".format(σ_st_st_left, μ_st_st_left, σ_st_st_right, μ_st_st_right))

                m_σ_left_range, m_σ_right_range = \
                    self._get_m_range_by_σμ_st_and_σμ_st_st(σ_st_st_left, isf_st_tuple_left[0], isf_st_tuple_right[0],
                                                            data_σ_μ.σ, data_sn, data_st)
                m_μ_left_range, m_μ_right_range = \
                    self._get_m_range_by_σμ_st_and_σμ_st_st(μ_st_st_left, isf_st_tuple_left[1], isf_st_tuple_right[1],
                                                            data_σ_μ.μ, data_sn, data_st)

                if _debug_condition(data_σ_μ):
                    logger.warning("m_σ_left_range, m_σ_right_range={}{}, | m_μ_left_range, m_μ_right_range={}{}"
                                   "".format(m_σ_left_range, m_σ_right_range, m_μ_left_range, m_μ_right_range))

                σ_μ_tn_element = 0
                σ_st_sn = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)]
                μ_st_sn = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)]

                if _debug_condition(data_σ_μ):
                    logger.warning("σ_st_sn={}, μ_st_sn={}".format(σ_st_sn, μ_st_sn))

                for m_σ_left, m_σ_right, m_μ_left, m_μ_right \
                        in product(range(*m_σ_left_range), range(*m_σ_right_range),
                                   range(*m_μ_left_range), range(*m_μ_right_range)):

                    if _debug_condition(data_σ_μ):
                        logger.warning("m_σ_left={}, m_σ_right={}, m_μ_left={}, m_μ_right={}"
                                       "".format(m_σ_left, m_σ_right, m_μ_left, m_μ_right))

                    σ_tn_element = σ_st_sn[m_σ_left - 1, m_σ_right - 1]
                    μ_tn_element = μ_st_sn[m_μ_left - 1, m_μ_right - 1]

                    if _debug_condition(data_σ_μ):
                        logger.warning("σ_tn_element={}, μ_tn_element={}"
                                       "".format(σ_tn_element, μ_tn_element))

                    σ_μ_tn_element += σ_tn_element * μ_tn_element

                # σ_tn_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_left - 1, m_σ_right - 1]
                # μ_tn_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_left - 1, m_μ_right - 1]

                if _debug_condition(data_σ_μ):
                    logger.warning("σ_μ_tn_element={}".format(σ_μ_tn_element))

                # σ_μ_tn_element = σ_tn_element * μ_tn_element
                # # τ看似没有拿，但实际已经包含在i, j里了
                # for left_j, right_j in product(range(len(isf_st_left_cols)), range(len(isf_st_right_cols))):
                #     isf_left_right = isf_st_left_square[left_i, left_j] * isf_st_right_square[right_i, right_j]
                #     matrix_element += σ_μ_tn_element * sp.sign(isf_left_right) * sp.sqrt(abs(isf_left_right))
                # isf_left_right = isf_st_left_square[left_i, right_i] * isf_st_right_square[right_i, left_i]
                # matrix_element += σ_μ_tn_element * sp.sign(isf_left_right) * sp.sqrt(abs(isf_left_right))
                isf_left_right = isf_st_left_square[left_i, left_i] * isf_st_right_square[right_i, right_i]
                matrix_element += σ_μ_tn_element * sp.sign(isf_left_right) * sp.sqrt(abs(isf_left_right))

        return True, matrix_element

    def _calc_isf_matrix_4_193d(self, row_index_list, ν_st, data_sn, data_σ_μ, data_st):
        """计算ISF的本征矩阵 选用4-193d优化
        4-193d的思想是用ISF代替CGC去构建矩阵"""
        if self.s_n <= 2:
            return self._calc_isf_matrix(row_index_list, ν_st, data_sn, data_σ_μ)

        matrix_div = len(row_index_list)
        isf_matrix = sp.zeros(matrix_div)
        h_ν_st = data_st.get_yt_num(ν_st)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("######## _calc_isf_matrix_4_193d")
            logger.warning("σ_({},{})={}".format(self.s_t, self.s_n, data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)]))
            logger.warning("μ_({},{})={}".format(self.s_t, self.s_n, data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)]))

        # ν_st的分支律
        ν_st_st_list = data_st.get_bl_yds(ν_st)
        for ν_st_st in ν_st_st_list:
            # m_ν_st_st = 1
            # m_ν_st = data_st.get_offset(ν_st_st, ν_st) + m_ν_st_st
            h_ν_st_st = data_st.get_yt_st_num(ν_st_st)
            sub_isf_matrix = sp.zeros(matrix_div)

            # 构建主对角元
            for i, row in enumerate(row_index_list):
                isf_st_tuple = (row[0], row[1], ν_st_st)

                if _debug_condition(data_σ_μ, ν_st):
                    logger.warning("isf_st_tuple={} for i={}".format(isf_st_tuple, i))

                flag, isf_matrix_element = \
                    self._calc_sub_isf_matrix_element_4_193d(isf_st_tuple, isf_st_tuple, data_sn, data_st, data_σ_μ)
                if not flag:
                    err_msg = "_calc_sub_isf_matrix_element_4_193d by isf_st_tuple={}, isf_st_tuple={}, " \
                              "data_sn={}, data_σ_μ={} with msg={}" \
                              "".format(isf_st_tuple, isf_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                sub_isf_matrix[i, i] = isf_matrix_element

            # 构建其他矩阵元素
            if matrix_div >= 2:
                for (left_i, left_st), (right_i, right_st) in combinations(enumerate(row_index_list), 2):
                    isf_st_tuple_left = (left_st[0], left_st[1], ν_st_st)
                    isf_st_tuple_right = (right_st[0], right_st[1], ν_st_st)

                    if _debug_condition(data_σ_μ, ν_st):
                        logger.warning("isf_st_tuple_left={}, isf_st_tuple_right={} for left_i={}, right_i={}"
                                       "".format(isf_st_tuple_left, isf_st_tuple_right, left_i, right_i))

                    flag, isf_matrix_element = \
                        self._calc_sub_isf_matrix_element_4_193d(isf_st_tuple_left, isf_st_tuple_right,
                                                                 data_sn, data_st, data_σ_μ)
                    if not flag:
                        err_msg = "_calc_sub_isf_matrix_element_4_193d by isf_st_tuple_left={}, " \
                                  "isf_st_tuple_right={}, data_sn={}, data_σ_μ={} with msg={}" \
                                  "".format(isf_st_tuple_left, isf_st_tuple_right, data_sn, data_σ_μ,
                                            isf_matrix_element)
                        logger.error(err_msg)
                        return False, err_msg
                    sub_isf_matrix[left_i, right_i] = isf_matrix_element
                    sub_isf_matrix[right_i, left_i] = isf_matrix_element  # 因为ISF的本征矩阵共轭

            # 融合
            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ ν_st_st={}, h_ν_st_st={}, sub_isf_matrix={}"
                               "".format(ν_st_st, h_ν_st_st, sub_isf_matrix))

            isf_matrix += h_ν_st_st * sub_isf_matrix

        isf_matrix = isf_matrix * (self.s_n - 1) / h_ν_st

        return True, isf_matrix

    def _calc_sym_isf_element(self, s_row, sym_σμν, sym_τ, sym_ν_st, sym_key, sym_d3, sym_k4,
                              data_sn, data_st, meta_data_σ_μ):
        """
        计算指定对称ISF element

        公式：
        这里，以meta_key=ν~σμ~为例 不带'是Sn的，带'是St的
        sqrt(h_ν'/h_ν) * ϵ(σ'μ'ν') * ϵ(σμν) * ISF_σσ'μμ'νν'
        = sqrt(h_μ'/h_μ) * ϵ(ν'~σ'μ'~) * ϵ(ν~σμ~) * Λ^ν_ν' * Λ^μ_μ' * ISF_ν~ν'~σσ'μ~μ'~
        有：
        ISF_ν~ν'~σσ'μ~μ'~
        = sqrt(h_ν'/h_ν)/sqrt(h_μ'/h_μ) * ϵ(σ'μ'ν')*ϵ(σμν) / ϵ(ν'~σ'μ'~)*ϵ(ν~σμ~) / Λ^ν_ν'*Λ^μ_μ' * ISF_σσ'μμ'νν'
        = sqrt(h_ν'*h_μ/h_ν*h_μ') * ϵ(σ'μ'ν')*ϵ(σμν) * ϵ(ν'~σ'μ'~)*ϵ(ν~σμ~) * Λ^ν_ν'*Λ^μ_μ' * ISF_σσ'μμ'νν'
        # 注意，因为实际输出的是ISF_square，所以sqrt是要平方的
        有：
        ISF_ν~ν'~σσ'μ~μ'~_square
        = h_ν'*h_μ / h_ν*h_μ' * ϵ(σ'μ'ν')*ϵ(σμν) * ϵ(ν'~σ'μ'~)*ϵ(ν~σμ~) * Λ^ν_ν'*Λ^μ_μ' * ISF_σσ'μμ'νν'_square
        """

        # if _debug_σμν_cond(*sym_σμν) and _debug_condition(meta_data_σ_μ):
        #     logger.warning("@@@@ catch sym_σμν={}, sym_ν_st={} with sym_τ={}, "
        #                    "s_row={}, sym_key={}, sym_d3={}, sym_k4={}"
        #                    "".format(sym_σμν, sym_ν_st, sym_τ, s_row, sym_key, sym_d3, sym_k4))

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
        ϵ_dict, _ = meta_data_σ_μ.get_ϵ_dict_and_flags(meta_σμν[-1], meta_τ)
        _, rst_st_dict = load_ϵ(*(self.s_t, *meta_σμν_st, meta_τ_st), is_with_flags=True)
        ϵ_st_dict, ϵ_st_flags = rst_st_dict.get("data", {}), rst_st_dict.get("flags", {})
        ϵ_4 = ϵ_st_dict["σμν"] * ϵ_dict["σμν"] * ϵ_st_dict[sym_key] * ϵ_dict[sym_key]
        # 3, Λ
        # 注意：这里的m是跟m'有关的，不是取自ϵ_flags
        # 法1，直接nb；法2，Λ^ν_ν' = Λ^ν(m) / Λ^ν'(m') = Λ^ν(m) * Λ^ν'(m')  # ：Λ^ν(m) = Λ^ν_ν' * Λ^ν'(m')
        meta_st_Λ_list_list = [data_st.get_phase_factor_list(yd) for yd in meta_σμν_st]
        meta_Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in meta_σμν]
        meta_m_st_list = ϵ_st_flags[sym_key]  # 其实取1也可以，m'是任意取的
        meta_m_list = [data_sn.quick_calc_m(m_st, yd_st, yd)
                       for m_st, yd_st, yd in zip(meta_m_st_list, meta_σμν_st, meta_σμν)]

        if _debug_σμν_cond(*sym_σμν) and _debug_condition(meta_data_σ_μ):
            logger.warning("@@@@ sym_key={}, ϵ_st={}, ϵ={}".format(sym_key, ϵ_st_dict[sym_key], ϵ_dict[sym_key]))
            logger.warning("@@@@ meta_m_st_list={}, meta_st_Λ_list_list={}".format(meta_m_st_list, meta_st_Λ_list_list))
            logger.warning("@@@@ meta_m_list={}, meta_Λ_list_list={}".format(meta_m_list, meta_Λ_list_list))

        ΛΛ = 1
        for d, k in zip(sym_d3, sym_k4):
            if k is True:
                ΛΛ *= meta_Λ_list_list[d][meta_m_list[d] - 1] * meta_st_Λ_list_list[d][meta_m_st_list[d] - 1]
        # 4, ISF_σσ'μμ'νν'
        meta_isf_square_dict = meta_data_σ_μ.get_isf_square_dict_by_ν_st(meta_σμν_st[-1])

        if _debug_σμν_cond(*sym_σμν) and _debug_condition(meta_data_σ_μ):
            logger.warning("@@@@ meta_σμν={}, meta_σμν_st={}, meta_isf_square_dict={}"
                           "".format(meta_σμν, meta_σμν_st, meta_isf_square_dict))
        meta_isf_rows = meta_isf_square_dict["rows"]  # rows有两种可能形式：1，([σ'], [μ'], τ')；2，([σ'], [μ'])
        meta_isf_cols = meta_isf_square_dict["cols"]  # cols有两种可能形式：1，[ν]；2，([ν], τ)
        meta_isf_square_matrix = meta_isf_square_dict["isf"]
        meta_row = (meta_σμν_st[0], meta_σμν_st[1]) if meta_τ_st is None \
            else (meta_σμν_st[0], meta_σμν_st[1], meta_τ_st)
        meta_col = meta_σμν[-1] if meta_τ is None else (meta_σμν[-1], meta_τ)

        if _debug_σμν_cond(*sym_σμν) and _debug_condition(meta_data_σ_μ):
            logger.warning("@@@@ meta_row={}, meta_isf_rows={}, meta_col={}, meta_isf_cols={}"
                           "".format(meta_row, meta_isf_rows, meta_col, meta_isf_cols))

        meta_row_index = meta_isf_rows.index(meta_row)
        meta_col_index = meta_isf_cols.index(meta_col)
        single_meta_isf_square = meta_isf_square_matrix[meta_row_index, meta_col_index]
        # ISF_ν~ν'~σσ'μ~μ'~
        single_sym_isf_square = h_square * ϵ_4 * ΛΛ * single_meta_isf_square

        if _debug_σμν_cond(*sym_σμν) and _debug_condition(meta_data_σ_μ):
            logger.warning("@@@@ single_sym_isf_square={} \n"
                           "= h_square={} * ϵ_4={} * ΛΛ={} * single_meta_isf_square={}"
                           "".format(single_sym_isf_square, h_square, ϵ_4, ΛΛ, single_meta_isf_square))

        return single_sym_isf_square

    def calc_σμ_symmetry_isf_include_save(self, meta_data_σ_μ, data_sn, data_st):
        """
        计算σμ全组的对称ISF
        注意：以ν结尾的ISF可以都算出来。而以σ、μ结尾的ISF只能算出部分列，这部分需要标记清楚

        1，收集当前σμ下所有ν的列表，为24个对称构型作准备（meta_data_σ_μ中已经有了，可以直接拿）
        2，记录24种ϵ下的所有新构型以及来源，all_json_sym_σμν_dict={json(σ_s,μ_s,ν_s): sym_key}
        3，根据all_json_sym_σμν_dict统计出不重不漏的σ_s,μ_s,ν_s'
        4，此时，ISF的title已经凑齐s_n,σ_s,μ_s,ν_s'，可以正式开启循环了
        5，首先计算rows，可以复用calc_row_indexes
        6，使用cg序列重建cols
        7，按列计算对称ISF的具体matrix
        8，保存
        """
        sym_isf_start_time = time.time()
        # 1，收集当前σμ下所有ν的列表，留下属于元σμν的部分，为24个对称构型作准备（meta_data_σ_μ中已经有了，可以直接拿）
        meta_ν_τ_list = meta_data_σ_μ.get_meta_ν_τ_list()
        # 2，读取元σμν经过24种ϵ变换后的所有新构型以及来源，all_json_sym_σμν_dict={(j_σ_s, j_μ_s, j_ν_s): sym_key}
        all_json_sym_σμν_dict = {}
        for meta_ν_τ in meta_ν_τ_list:
            meta_ν, meta_τ = meta_ν_τ if isinstance(meta_ν_τ, tuple) else (meta_ν_τ, None)
            if meta_τ is not None and meta_τ != 1:  # 只留None和1
                continue
            meta_σμν = (meta_data_σ_μ.σ, meta_data_σ_μ.μ, meta_ν)
            _, json_sym_σμν_dict = load_sym_σμν(*(self.s_n, *meta_σμν), is_flag_true_if_not_s_n=False)
            for j_sym_σμν, sym_key_list in json_sym_σμν_dict.items():
                if j_sym_σμν not in all_json_sym_σμν_dict:
                    all_json_sym_σμν_dict[j_sym_σμν] = sym_key_list[0]  # 记录一个就够了
        # 3，根据all_json_sym_σμν_dict统计出不重不漏的σ_s,μ_s,ν_s'
        json_sym_σμ_ν_st_dict = {}  # {(j_σ_s, j_μ_s): [ν_s']}
        for (j_sym_σ, j_sym_μ, j_sym_ν) in all_json_sym_σμν_dict.keys():
            if (j_sym_σ, j_sym_μ) not in json_sym_σμ_ν_st_dict:
                json_sym_σμ_ν_st_dict[(j_sym_σ, j_sym_μ)] = []
            sym_ν_bl = data_sn.get_bl_yds(json.loads(j_sym_ν))
            for sym_ν_st in sym_ν_bl:
                if sym_ν_st not in json_sym_σμ_ν_st_dict[(j_sym_σ, j_sym_μ)]:
                    json_sym_σμ_ν_st_dict[(j_sym_σ, j_sym_μ)].append(sym_ν_st)

        if _debug_condition(meta_data_σ_μ):
            logger.warning("@@@@ meta_data_σ_μ={}, meta_ν_τ_list={}, "
                           "all_json_sym_σμν_dict={}, json_sym_σμ_ν_st_dict={}"
                           "".format(meta_data_σ_μ, meta_ν_τ_list, all_json_sym_σμν_dict, json_sym_σμ_ν_st_dict))

        # 4，此时，ISF的title已经凑齐s_n,σ_s,μ_s,ν_s'，可以正式开启循环了
        sym_data_σ_μ_flag = (None, None)
        sym_data_σ_μ = None
        for (j_sym_σ, j_sym_μ), sym_ν_st_list in json_sym_σμ_ν_st_dict.items():
            sym_σ, sym_μ = map(json.loads, [j_sym_σ, j_sym_μ])
            for sym_ν_st in sym_ν_st_list:
                # 正式拿到s_n,σ_s,μ_s,ν_s'
                sym_isf_tuple = (self.s_n, sym_σ, sym_μ, sym_ν_st)
                _, old_isf_info_dict = load_isf(*sym_isf_tuple, output_mode="all", ex_params=["data", "flags"],
                                                is_flag_true_if_not_s_n=True)
                if old_isf_info_dict and old_isf_info_dict.get("flags", {}).get("finish_cols", None) is None:
                    continue
                if sym_data_σ_μ_flag != (sym_σ, sym_μ):
                    sym_data_σ_μ = SimpleΣMDataHelper(self.s_n, sym_σ, sym_μ, data_sn, data_st)
                    sym_data_σ_μ_flag = (sym_σ, sym_μ)
                if old_isf_info_dict is not False:
                    old_isf_info_dict_cp = copy.deepcopy(old_isf_info_dict)
                    sym_isf_square_dict = old_isf_info_dict_cp["data"]
                    sym_rows = sym_isf_square_dict["rows"]
                    sym_cols = sym_isf_square_dict["cols"]
                    sym_isf_square_matrix = sym_isf_square_dict["isf"]
                    finish_cols = old_isf_info_dict_cp["flags"]["finish_cols"]
                    old_adding_time = old_isf_info_dict_cp["flags"]["speed_time"]
                else:
                    # 5，首先计算rows，可以复用calc_row_indexes
                    # 注意，由meta_rows直接变换的话，值相同，但顺序不同！
                    sym_rows = self.calc_row_indexes(sym_data_σ_μ, sym_ν_st, data_st)
                    # 6，使用cg序列重建cols
                    sym_cols = self.calc_col_indexes(sym_data_σ_μ, sym_ν_st, data_sn)
                    # 7，按列计算对称ISF的具体matrix
                    if len(sym_rows) != len(sym_cols):
                        err_msg = "len(sym_rows)={} len(sym_cols)={} should eq but not " \
                                  "with meta_data_σ_μ={}, sym_data_σ_μ={}" \
                                  "".format(len(sym_rows), len(sym_cols), meta_data_σ_μ, sym_data_σ_μ)
                        return False, err_msg
                    sym_div = len(sym_rows)
                    if sym_div == 0:
                        continue
                    sym_isf_square_matrix = sp.zeros(sym_div)
                    finish_cols = []
                    old_adding_time = 0
                    sym_isf_square_dict = {"rows": sym_rows,
                                           "cols": sym_cols,
                                           "isf": sym_isf_square_matrix}

                    if _debug_condition(meta_data_σ_μ) and _debug_isf_cond(sym_σ, sym_μ, sym_ν_st):
                        logger.warning("@@@@ sym_σ={}, sym_μ={}, sym_ν_st={}, "
                                       "old_isf_info_dict={}, sym_isf_square_dict={}"
                                       "".format(sym_σ, sym_μ, sym_ν_st,
                                                 old_isf_info_dict, sym_isf_square_dict))

                for s_col_i, s_col in enumerate(sym_cols):
                    if s_col in finish_cols:
                        continue
                    sym_ν, sym_τ = s_col if isinstance(s_col, tuple) else (s_col, None)
                    sym_σμν = (sym_σ, sym_μ, sym_ν)
                    json_sym_σμν_key = (j_sym_σ, j_sym_μ, json.dumps(sym_ν))
                    if json_sym_σμν_key not in all_json_sym_σμν_dict:
                        continue  # 当前的元σμν一般不能通过对称得到当前所有isf的ν（不能的那部分，就是当前还未循环到的σμ组合及其对称）
                    finish_cols.append(s_col)
                    sym_key = all_json_sym_σμν_dict[json_sym_σμν_key]
                    sym_d3, sym_k4 = data_sn.get_d3_k4_by_ϵ_key(sym_key)

                    if _debug_condition(meta_data_σ_μ) and _debug_isf_cond(sym_σ, sym_μ, sym_ν_st):
                        logger.warning("@@@@ before calc s_col={}, sym_isf_square_matrix={}"
                                       "".format(s_col, sym_isf_square_matrix))

                    for s_row_i, s_row in enumerate(sym_rows):
                        sym_isf_element = self._calc_sym_isf_element(s_row, sym_σμν, sym_τ, sym_ν_st, sym_key,
                                                                     sym_d3, sym_k4, data_sn, data_st, meta_data_σ_μ)
                        sym_isf_square_matrix[s_row_i, s_col_i] = sym_isf_element

                        if _debug_condition(meta_data_σ_μ) and _debug_isf_cond(sym_σ, sym_μ, sym_ν_st) \
                                and s_col == [2, 1, 1]:
                            logger.warning("@@@@ test 'μν~σ~' for _calc_sym_isf_element again but not save")
                            test_isf_element = self._calc_sym_isf_element(s_row, sym_σμν, sym_τ, sym_ν_st, 'μν~σ~',
                                                                          group_d3[5], group_k4[3], data_sn, data_st,
                                                                          meta_data_σ_μ)
                            logger.warning("@@@@ sym_isf_element={}, test_isf_element={}"
                                           "".format(sym_isf_element, test_isf_element))

                    if _debug_condition(meta_data_σ_μ) and _debug_isf_cond(sym_σ, sym_μ, sym_ν_st):
                        logger.warning("@@@@ after calc s_col={}, sym_isf_square_matrix={}"
                                       "".format(s_col, sym_isf_square_matrix))

                # 8，保存
                # 根据python特性，sym_isf_square_dict是浅拷贝，对其元素所新指定的变量赋值，就相当于更新了它

                if _debug_condition(meta_data_σ_μ) and _debug_isf_cond(sym_σ, sym_μ, sym_ν_st):
                    logger.warning("@@@@ meta_data_σ_μ={}, sym_σ={}, sym_μ={}, sym_ν_st={}, "
                                   "sym_isf_square_dict={}, finish_cols={}"
                                   "".format(meta_data_σ_μ, sym_σ, sym_μ, sym_ν_st, sym_isf_square_dict, finish_cols))

                sym_isf_speed_time = int(time.time() - sym_isf_start_time)
                if len(finish_cols) == len(sym_cols):
                    flag, msg = save_isf(self.s_n, sym_σ, sym_μ, sym_ν_st, sym_isf_square_dict,
                                         sym_isf_speed_time + old_adding_time)
                    if not flag:
                        err_msg = "save meta isf meet error with s_i={}, sym_σ={}, sym_μ={}, sym_ν_st={}, " \
                                  "sym_isf_square_dict={}, " \
                                  "msg={}".format(self.s_n, sym_σ, sym_μ, sym_ν_st, sym_isf_square_dict, msg)
                        logger.error(err_msg)
                        return False, err_msg
                else:
                    flag, msg = save_isf(self.s_n, sym_σ, sym_μ, sym_ν_st, sym_isf_square_dict,
                                         sym_isf_speed_time + old_adding_time, finish_cols=finish_cols)
                    if not flag:
                        err_msg = "save meta isf meet error with s_i={}, sym_σ={}, sym_μ={}, sym_ν_st={}, " \
                                  "sym_isf_square_dict={}, finish_cols={}, msg={}" \
                                  "".format(self.s_n, sym_σ, sym_μ, sym_ν_st, sym_isf_square_dict, finish_cols, msg)
                        logger.error(err_msg)
                        return False, err_msg
                sym_isf_start_time = time.time()

        return True, None

    # @staticmethod
    # def _sort_isf_rows(symmetry_isf_rows_unsort, data_st):
    #     """将未按照Yamanouchi排序的symmetry_isf_rows_unsort按照Yamanouchi排序
    #     技巧是利用yd_st_index达到排序yd_st的目的"""
    #     sort_isf_rows = []
    #     sym_2_meta_relationship = []  # 形如[1, 3, 2]表示已排序的rows是从meta中哪里来的
    #     all_yd_st_list = data_st.yd_list
    #     σ_st_index_list = [all_yd_st_list.index(i[0]) for i in symmetry_isf_rows_unsort]
    #     σ_st_index_one_list = list(set(σ_st_index_list))
    #     for σ_st_index_one in sorted(σ_st_index_one_list):  # 排序金牌的σ
    #         σ_st = all_yd_st_list[σ_st_index_one]
    #         μ_st_index_list = [all_yd_st_list.index(i[1]) for i in symmetry_isf_rows_unsort if i[0] == σ_st]
    #         μ_st_index_one_list = list(set(μ_st_index_list))
    #         for μ_st_index_one in sorted(μ_st_index_one_list):  # 排序银牌的μ
    #             μ_st = all_yd_st_list[μ_st_index_one]
    #             # 实有τ的才会是real list，没有的会是[]
    #             τ_st_list = [i[2] for i in symmetry_isf_rows_unsort if len(i) == 3 and i[0] == σ_st and i[1] == μ_st]
    #             if τ_st_list:
    #                 for τ_st in sorted(τ_st_list):  # 排序铜牌的τ
    #                     single_tuple = (σ_st, μ_st, τ_st)
    #                     sort_isf_rows.append(single_tuple)
    #                     sym_2_meta_relationship.append(symmetry_isf_rows_unsort.index(single_tuple))
    #             else:
    #                 single_tuple = (σ_st, μ_st)
    #                 sort_isf_rows.append(single_tuple)
    #                 sym_2_meta_relationship.append(symmetry_isf_rows_unsort.index(single_tuple))
    #     # 形如[1, 3, 2]表示meta会到已排序的rows哪里去
    #     meta_2_sym_relationship = [sym_2_meta_relationship.index(i) for i in range(len(sym_2_meta_relationship))]
    #     return sort_isf_rows, meta_2_sym_relationship, sym_2_meta_relationship

    # @staticmethod
    # def _sort_isf_cols(tilde_isf_cols_unsort, data_sn):
    #     """将未按照Yamanouchi排序的tilde_isf_cols_unsort按照Yamanouchi排序
    #     技巧是利用yd_index达到排序yd的目的"""
    #     sort_isf_cols = []
    #     sym_2_meta_relationship = []  # 形如[1, 3, 2]表示已排序的cols是从meta中哪里来的
    #     all_yd_list = data_sn.yd_list
    #     ν_index_list = [all_yd_list.index(i[0]) if isinstance(i, tuple) else all_yd_list.index(i)
    #                     for i in tilde_isf_cols_unsort]
    #     ν_index_one_list = list(set(ν_index_list))
    #     for ν_index_one in sorted(ν_index_one_list):  # 排序金牌的ν
    #         ν = all_yd_list[ν_index_one]
    #         # 实有τ的才会是real list，没有的会是[]
    #         τ_list = [i[1] for i in tilde_isf_cols_unsort if isinstance(i, tuple) and i[0] == ν]
    #         if τ_list:
    #             for τ in sorted(τ_list):  # 排序银牌的τ
    #                 single_tuple = (ν, τ)
    #                 sort_isf_cols.append(single_tuple)
    #                 sym_2_meta_relationship.append(tilde_isf_cols_unsort.index(single_tuple))
    #         else:
    #             sort_isf_cols.append(ν)
    #             sym_2_meta_relationship.append(tilde_isf_cols_unsort.index(ν))
    #     # 形如[1, 3, 2]表示meta会到已排序的cols哪里去
    #     meta_2_sym_relationship = [sym_2_meta_relationship.index(i) for i in range(len(sym_2_meta_relationship))]
    #     return sort_isf_cols, meta_2_sym_relationship, sym_2_meta_relationship


class CGCHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(CGCHelper, self).__init__()

    def calc_cgc_by_isf_include_save(self, data_σ_μ, ν_st, data_sn, data_st, non_meta_ν_end=None):
        """
        使用ISF计算CGC
        
        公式：4-189a、4-195a
        
        当new_ν=None时：
        凭借ISF中的元那部分，计算和储存元CGC

        当non_meta_ν_end=yd时：
        利用ISF计算和储存σνμ和νμσ的CGC
        """
        # 准备数据
        isf_square_dict = data_σ_μ.get_isf_square_dict_by_ν_st(ν_st)
        isf_rows = isf_square_dict["rows"]
        isf_cols = isf_square_dict["cols"]
        isf_square_matrix = isf_square_dict["isf"]
        isf_finish_cols = data_σ_μ.get_isf_finish_cols_by_ν_st(ν_st)
        isf_finish_ν_list = []
        for col in isf_finish_cols:
            ν, τ = col if isinstance(col, tuple) else (col, None)
            if ν not in isf_finish_ν_list:
                isf_finish_ν_list.append(ν)
        # 算metaCGC的时候，只留meta_ν_list_of_σμ与isf_finish_ν_list重合那部分col；算non_meta_ν_end的时候，只留等于它的col
        if non_meta_ν_end is not None:
            wanted_ν_list = [non_meta_ν_end]
        else:
            wanted_ν_list = []
            for meta_ν in data_σ_μ.meta_ν_list_of_σμ:
                if meta_ν in isf_finish_ν_list:
                    wanted_ν_list.append(meta_ν)
        # 得到h_ν_st
        h_ν_st = data_st.get_yt_num(ν_st)

        # catch_cond_31 = (data_σ_μ.σ == data_σ_μ.μ == [3, 1])
        # if catch_cond_31:
        #     logger.warning("######## catch CGC maker [3, 1] * [3, 1], "
        #                    "isf_finish_ν_list={}, wanted_ν_list={}, h_ν_st={}"
        #                    "".format(isf_finish_ν_list, wanted_ν_list, h_ν_st))

        # 下面的循环，由互不相干的四层构成：
        # m_ν'（CGC名）；ν,τ（ISF列）；ISF行、(m_σ', m_μ')(CGC内)。更改它们的顺序只影响算法速度，不影响计算结果。
        for m_ν_st in range(1, h_ν_st + 1):  # m_ν'循环
            # 在深层循环前，准备数据，空间换时间
            cgc_st_square_dict_list_by_row = []
            offset_σ_list = []
            offset_μ_list = []
            for isf_row in isf_rows:  # 数据准备，不是四层之一
                σ_st, μ_st, τ_st = isf_row if len(isf_row) == 3 else (*isf_row, None)
                flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st,
                                                    is_flag_true_if_not_s_n=False)
                if not flag:
                    err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, τ_st={}, " \
                              "m_ν_st={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st,
                                                         cgc_st_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
                cgc_st_square_dict_list_by_row.append(cgc_st_square_dict)
                offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
                offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
                offset_σ_list.append(offset_σ)
                offset_μ_list.append(offset_μ)

            # if catch_cond_31:
            #     logger.warning("## isf_rows={}, m_ν_st={}, ν_st={}, "
            #                    "cgc_st_square_dict_list_by_row={}, offset_σ_list={}, offset_μ_list={}"
            #                    "".format(isf_rows, m_ν_st, ν_st,
            #                              cgc_st_square_dict_list_by_row, offset_σ_list, offset_μ_list))

            for col_i, isf_col in enumerate(isf_cols):  # 按照列循环
                ν, τ = isf_col if isinstance(isf_col, tuple) else (isf_col, None)

                # if catch_cond_31:
                #     logger.warning("#1# ν={}, τ={}".format(ν, τ))

                if ν not in wanted_ν_list:
                    continue

                # 1.2，根据元ISF计算元CGC【[σ][μ][ν]τ】
                # 拆解列信息
                isf_square_vector = isf_square_matrix.col(col_i)
                # 每个ν有自己的offset
                offset_ν = data_sn.get_offset(ν_st, ν)
                m_ν = offset_ν + m_ν_st
                # 1.2.1，已经得到σ,μ,ν,τ,m，开始计算元CGC前先检查一下（因为非元σμν那部分已经算好了）
                _, is_calc_ed = is_cgc_exist(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν)
                if is_calc_ed is True:
                    continue
                # 1.2.2，判断，并注册元信息
                data_σ_μ.register_meta_ν_τ_list(ν, τ)
                meta_cgc_start_time = time.time()
                meta_cgc_square_dict = {}

                # if catch_cond_31:
                #     logger.warning("#2# ν={}, τ={}".format(ν, τ))
                #
                # catch_cond_31_22 = (catch_cond_31 and ν == [2, 2])
                #
                # if catch_cond_31:
                #     logger.warning("#3# ν={}, τ={}".format(ν, τ))
                #
                # if catch_cond_31_22:
                #     logger.warning("#4# ν={}, τ={}".format(ν, τ))
                #
                # if catch_cond_31_22:
                #     logger.warning("## catch CGC maker [3, 1] * [3, 1] = [2, 2] "
                #                    "isf_col={}, ν={}, τ={}, m_ν_st={}"
                #                    "".format(isf_col, ν, τ, m_ν_st))

                for isf_row, isf_square, single_cgc_st_square_dict, single_offset_σ, single_offset_μ in \
                        zip(isf_rows, isf_square_vector, cgc_st_square_dict_list_by_row,
                            offset_σ_list, offset_μ_list):  # 按照行循环
                    if isf_square == 0:
                        continue

                    # if catch_cond_31_22:
                    #     logger.warning("## isf_row={}, isf_square={}, single_cgc_st_square_dict={}, "
                    #                    "single_offset_σ={}, single_offset_μ={}"
                    #                    "".format(isf_row, isf_square, single_cgc_st_square_dict,
                    #                              single_offset_σ, single_offset_μ))

                    for (m_σ_st, m_μ_st), cgc_st_square in single_cgc_st_square_dict.items():  # (m_σ', m_μ')循环
                        if cgc_st_square == 0:
                            continue
                        cgc_key = (single_offset_σ + m_σ_st, single_offset_μ + m_μ_st)

                        # if catch_cond_31_22:
                        #     logger.warning("## m_σ_st={}, m_μ_st={}, cgc_st_square={}, cgc_key={}"
                        #                    "".format(m_σ_st, m_μ_st, cgc_st_square, cgc_key))

                        if cgc_key not in meta_cgc_square_dict:
                            single_cgc_square = isf_square * cgc_st_square
                            if single_cgc_square != 0:
                                meta_cgc_square_dict[cgc_key] = isf_square * cgc_st_square
                        else:
                            cgc_new_part = sp.sign(isf_square) * sp.sign(cgc_st_square) \
                                           * sp.sqrt(abs(isf_square * cgc_st_square))
                            cgc_old_part = sp.sign(meta_cgc_square_dict[cgc_key]) \
                                           * sp.sqrt(abs(meta_cgc_square_dict[cgc_key]))
                            update_cgc_square = sp.sign(cgc_new_part + cgc_old_part) \
                                                * (cgc_new_part + cgc_old_part)**2
                            if update_cgc_square != 0:
                                meta_cgc_square_dict[cgc_key] = update_cgc_square  # 覆盖
                            else:
                                meta_cgc_square_dict.pop(cgc_key)

                        if self.s_n == 6 and data_σ_μ.σ == [3, 2, 1] and data_σ_μ.μ == [3, 2, 1] \
                                and ν == [3, 2, 1] and τ == 1 and m_ν == 1 and cgc_key == (1, 1):
                            logger.warning("$$$$ [3, 2, 1]*[3, 2, 1]=[3, 2, 1]τ5m1: (1, 3): {} from\n"
                                           "isf: ν_st={}, isf_row={}, isf_col={}, isf_square={}\n"
                                           "cgc: m_σ_st={}, m_μ_st={}, single_cgc_st_square_dict={}"
                                           "".format(meta_cgc_square_dict[cgc_key],
                                                     ν_st, isf_row, isf_col, isf_square,
                                                     m_σ_st, m_μ_st, single_cgc_st_square_dict))

                # if catch_cond_31_22:
                #     logger.warning("## meta_cgc_square_dict={}, m_ν={}".format(meta_cgc_square_dict, m_ν))

                cgc_square_n = sum(abs(i) for i in meta_cgc_square_dict.values())
                if cgc_square_n != 1:  # 归一性检查  # 由4-189a，循环完所有τ'，CGC理论上就是完成了的
                    err_msg = "calc cgc fail by self.s_n={}, σ={}, μ={}, " \
                              "ν={}, τ={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "meet cgc_square_n={} not eq 1".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν,
                                                                     meta_cgc_square_dict, cgc_square_n)
                    logger.error(err_msg)
                    return False, err_msg

                meta_cgc_speed_time = int(time.time() - meta_cgc_start_time)
                flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, meta_cgc_square_dict,
                                     meta_cgc_speed_time)
                if not flag:
                    err_msg = "save_cgc fail by self.s_n={}, σ={}, μ={}, " \
                              "ν={}, τ={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, meta_cgc_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg
                # 对于non_meta_ν_end=yd的情况，在这里继续下去，就不用再以后重新load一次CGC了，但为了复用代码，我们浪费一点点资源

        return True, None

    def calc_sym_cgc_by_same_end_cgc_include_save(self, current_data_σ_μ, current_ν, τ, data_sn):
        """
        利用一组当前已知CGC和同它结尾相同的ϵ，计算同它结尾相同的对称CGC

        公式：
        CGC^sym_mode(m|meta_key cond)
        = ϵ(σμν) * h_σ/h_ν * ϵ(sym_mode) * ΛΛ_u^sym_mode(m｜meta_key cond) * CGC^meta_mode[m｜meta_key cond]
        因为只考虑ν结尾的情况，所有，有：
        = ϵ(σμν) * ϵ(sym_mode) * ΛΛ_u^sym_mode(m｜meta_key cond) * CGC^meta_mode[m｜meta_key cond]
        """
        # 数据准备
        sym_cgc_start_time = time.time()
        current_yd_σμν = (current_data_σ_μ.σ, current_data_σ_μ.μ, current_ν)
        current_Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in current_yd_σμν]
        current_h_list = [data_sn.get_yt_num(yd) for yd in current_yd_σμν]
        ϵ_yds_dict = {}  # {ϵ_key: sym_current_yd_σμν}
        for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
            sym_current_yd_σμν = tuple(current_yd_σμν[d] if k is False else
                                       data_sn.get_tilde(current_yd_σμν[d]) for d, k in zip(d3, k4))
            ϵ_yds_dict[ϵ_key] = sym_current_yd_σμν
        # 对于meta，current_ϵ_dict必然能拿到；对于非meta，在使用前，有代码负责提前装载，所以也能拿到
        # current_ϵ_dict, current_ϵ_flags = current_data_σ_μ.get_ϵ_dict_and_flags(current_ν, τ)
        current_ϵ_dict, _ = current_data_σ_μ.get_ϵ_dict_and_flags(current_ν, τ)
        # ϵ(σμν) * ϵ(sym_mode) 与m无关，可以拿到前面算
        ϵ_2_dict = {}
        for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
            ϵ_2 = current_ϵ_dict["σμν"] * current_ϵ_dict[ϵ_key]
            ϵ_2_dict[ϵ_key] = ϵ_2
        h_ν = data_sn.get_yt_num(current_ν)
        for current_m_ν in range(1, h_ν + 1):
            current_m_ν_tilde = h_ν + 1 - current_m_ν
            # 准备结果的空字典  # 对于一个m，一次循环里，把所有ν结尾的对称CGC都算出来（最多也就7个），空间换时间
            sym_cgc_square_dict_dict = {}  # {ϵ_key: sym_cgc_square_dict}
            sym_cgc_tuple_list = []  # 记录要算哪些tuple，避免重复
            for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                sym_current_yd_σμν = ϵ_yds_dict[ϵ_key]
                cgc_tuple = (self.s_n, *sym_current_yd_σμν, τ, current_m_ν) if k4[-1] is False \
                    else (self.s_n, *sym_current_yd_σμν, τ, current_m_ν_tilde)
                _, is_calc_ed = is_cgc_exist(*cgc_tuple)
                if is_calc_ed is False and cgc_tuple not in sym_cgc_tuple_list:  # 这样，就可以保证不重不漏了
                    sym_cgc_tuple_list.append(cgc_tuple)
                    sym_cgc_square_dict_dict[ϵ_key] = {}
            sym_cgc_square_dict_dict_keys = list(sym_cgc_square_dict_dict.keys())
            if not sym_cgc_square_dict_dict:  # 空字典说明对称都被算过了
                continue

            # catch_cond = (4, [2, 1, 1], [3, 1], [2, 2], None, 1) in sym_cgc_tuple_list or \
            #              (4, [2, 1, 1], [3, 1], [2, 2], None, 2) in sym_cgc_tuple_list
            # catch_cond = (3, [2, 1], [3], [2, 1], None, 1) in sym_cgc_tuple_list or \
            #              (3, [2, 1], [3], [2, 1], None, 2) in sym_cgc_tuple_list
            # if catch_cond:
            #     logger.warning("######## catch (3, [2, 1], [3], [2, 1]) with sym_cgc_square_dict_dict={}, "
            #                    "current_yd_σμν={}, τ={}, current_m_ν={}, sym_cgc_tuple_list={}"
            #                    "".format(sym_cgc_square_dict_dict, current_yd_σμν, τ, current_m_ν, sym_cgc_tuple_list))

            # 准备meta_cgc
            flag, meta_cgc_square_dict = load_cgc(*(self.s_n, *current_yd_σμν, τ, current_m_ν),
                                                  is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_cgc fail by self.s_n={}, current_yd_σμν={}, τ={}, current_m_ν={}, " \
                          "with msg={}".format(self.s_n, current_yd_σμν, τ, current_m_ν, meta_cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg
            '''这里，meta_cgc_square_dict是最大的，所以它的循环应该放在外面，循环最少的次数，而且它循环，可以在外面就尽量在外面'''
            for (meta_m_σ, meta_m_μ), meta_cgc_square in meta_cgc_square_dict.items():
                meta_m_σμν = (meta_m_σ, meta_m_μ, current_m_ν)
                for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                    if ϵ_key not in sym_cgc_square_dict_dict_keys:
                        continue
                    # ϵ(σμν) * ϵ(sym_mode)已经在外面算好了
                    ϵ_2 = ϵ_2_dict[ϵ_key]
                    # ΛΛ_u^sym_mode(m｜meta_key cond) & sym m
                    ΛΛ = 1
                    sym_m_list = [meta_m_σμν[d] if k is False else
                                  current_h_list[d] + 1 - meta_m_σμν[d] for d, k in zip(d3, k4)]
                    for i, d, k in zip(range(3), d3, k4):
                        if k is True:
                            ΛΛ *= current_Λ_list_list[d][meta_m_σμν[d] - 1]
                    # 小心，下面注释掉的是错误代码，它错在拿了meta_m_σ，却去算sym_m_σ
                    # for d, k, single_meta_m_σμν, single_meta_h in zip(d3, k4, meta_m_σμν, current_h_list):
                    #     if k is True:
                    #         ΛΛ *= current_Λ_list_list[d][meta_m_σμν[d] - 1]
                    #         sym_m_list[d] = single_meta_h + 1 - single_meta_m_σμν
                    #     else:
                    #         sym_m_list[d] = single_meta_m_σμν

                        # if catch_cond and (ϵ_key == "μσν" or ϵ_key == "μ~σν~"):
                        #     logger.warning("ϵ_key={}, d={}, k={}, sym_m_list[i]={}"
                        #                    "".format(ϵ_key, d, k, sym_m_list[i]))

                    # CGC^meta_mode[m｜meta_key cond]
                    # meta_cgc_square
                    # CGC^sym_mode(m|meta_key cond)
                    sym_cgc_square = ϵ_2 * ΛΛ * meta_cgc_square
                    sym_cgc_square_dict_dict[ϵ_key][(sym_m_list[0], sym_m_list[1])] = sym_cgc_square

                    # if catch_cond and ϵ_key == "σ~μν~":
                    # if catch_cond and (ϵ_key == "μσν" or ϵ_key == "μ~σν~"):
                    #     logger.warning("## ϵ_key={}, sym_cgc_square={} eq\nϵ_2={} * ΛΛ={} * meta_cgc_square={} "
                    #                    "with sym_m_list={}, meta_m_σμν={}, sym_cgc_square_dict_dict={}"
                    #                    "".format(ϵ_key, sym_cgc_square, ϵ_2, ΛΛ, meta_cgc_square, sym_m_list,
                    #                              meta_m_σμν, sym_cgc_square_dict_dict))

            # if catch_cond:
            #     logger.warning("## sym_cgc_square_dict_dict={}, sym_cgc_tuple_list={}, current_m_ν={}"
            #                    "".format(sym_cgc_square_dict_dict, sym_cgc_tuple_list, current_m_ν))

            # 保存  # 这里的时间采用平均值来计算
            sym_cgc_speed_time = int(time.time() - sym_cgc_start_time)
            sym_cgc_speed_time_average = int(sym_cgc_speed_time / len(sym_cgc_square_dict_dict_keys))
            sym_cgc_start_time = time.time()
            for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                if ϵ_key not in sym_cgc_square_dict_dict_keys:
                    continue
                sym_current_yd_σμν = ϵ_yds_dict[ϵ_key]
                cgc_tuple = (self.s_n, *sym_current_yd_σμν, τ, current_m_ν) if k4[-1] is False \
                    else (self.s_n, *sym_current_yd_σμν, τ, current_m_ν_tilde)
                flag, msg = save_cgc(*(*cgc_tuple, sym_cgc_square_dict_dict[ϵ_key], sym_cgc_speed_time_average))
                if not flag:
                    err_msg = "save_cgc fail by cgc_tuple={}, sym_cgc_square_dict_dict[ϵ_key]={}, ϵ_key={} with " \
                              "msg={}".format(cgc_tuple, sym_cgc_square_dict_dict[ϵ_key], ϵ_key, msg)
                    logger.error(err_msg)
                    return False, err_msg

        return True, None

    def calc_sym_μ_σ_end_cgc_include_save(self, new_σμν, data_sn, data_st):
        """
        1，利用ISF计算σνμ或νμσ的CGC
        2，以σνμ或νμσ为起点，对称出同样以μ或σ结尾的CGC
        """
        # 准备
        new_σ, new_μ, new_ν = new_σμν
        new_data_σ_μ = ΣMDataHelper(self.s_n, new_σ, new_μ, data_sn, data_st)
        # new_data_σ_μ.hook_meta_ν_list_of_σμ_and_get_ν_st_list_inner_meta_bl()  # 这里按定义没有meta可挂载
        cg_series_list = new_data_σ_μ.cg_series_list
        yd_list = data_sn.yd_list
        new_τ_max = cg_series_list[yd_list.index(new_ν)]  # 预先拿到new_τ_max
        for new_τ in range(1, new_τ_max + 1):
            new_τ = None if new_τ_max == 1 else new_τ
            _, new_ϵ = load_ϵ(*(self.s_n, *new_σμν, new_τ), is_with_flags=True)
            new_ϵ_dict = new_ϵ["data"]
            new_ϵ_flags = new_ϵ["flags"]
            new_data_σ_μ.register_ϵ_dict_and_flags(new_ν, new_τ, new_ϵ_dict, new_ϵ_flags)  # 预先装载ϵ
        # new_meta_ν的所有St分支，必然包含new_meta_ν的ISF，也必然已经被计算了（指的是ISF中含有new_meta_ν的列）
        for new_ν_st in data_sn.get_bl_yds(new_ν):
            _, is_calc_ed = is_isf_exist(*(self.s_n, new_σ, new_μ, new_ν_st))
            if is_calc_ed is False:
                err_msg = "isf must exist but not with s_n={}, new_σμν={}, new_ν_st={}".format(
                    self.s_n, new_σμν, new_ν_st)
                logger.error(err_msg)
                return False, err_msg
            # 1，利用ISF计算σνμ或νμσ的CGC
            # 对于partial_isf_square_dict，我们只取用new_ν，有就行。其他部分是否完整并不关心
            _, new_ν_isf_info_dict = load_isf(self.s_n, new_σ, new_μ, new_ν_st,
                                              output_mode="all", ex_params=["data", "flags"])
            if new_ν_isf_info_dict.get("flags", {}).get("finish_cols", None) is None:
                new_ν_finish_cols = new_ν_isf_info_dict["data"]["cols"]
            else:
                new_ν_finish_cols = new_ν_isf_info_dict["flags"]["finish_cols"]
            new_data_σ_μ.register_isf_square_dict_and_finish_cols_by_ν_st(new_ν_st, new_ν_isf_info_dict["data"],
                                                                          new_ν_finish_cols)

            # 这里，我们以多load一次new_σμν所有CGC的代价，复用calc_meta_cgc_include_save
            flag, msg = self.calc_cgc_by_isf_include_save(new_data_σ_μ, new_ν_st, data_sn, data_st,
                                                          non_meta_ν_end=new_ν)
            if not flag:
                err_msg = "calc_cgc_by_isf_include_save fail by new_data_σ_μ={}, new_ν_st={}, data_sn={}, data_st={}" \
                          ", non_meta_ν_end={}, with msg={}" \
                          "".format(new_data_σ_μ, new_ν_st, data_sn, data_st, new_ν, msg)
                logger.error(err_msg)
                return False, err_msg

        # 2，以σνμ或νμσ为起点，对称出同样以μ或σ结尾的CGC
        for new_τ in range(1, new_τ_max + 1):
            new_τ = None if new_τ_max == 1 else new_τ
            flag, msg = self.calc_sym_cgc_by_same_end_cgc_include_save(new_data_σ_μ, new_ν, new_τ, data_sn)
            if not flag:
                err_msg = "calc_sym_cgc_by_same_end_cgc_include_save fail by new_data_σ_μ={}, new_ν={}, " \
                          "new_τ={}, data_sn={} with msg={}".format(new_data_σ_μ, new_ν, new_τ, data_sn, msg)
                logger.error(err_msg)
                return False, err_msg

        return True, None


class EHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(EHelper, self).__init__()

    def calc_meta_ϵ_dict(self, data_σ_μ, ν, τ, data_sn):
        """根据定义计算6组24个ϵ
        ϵ的意义最重要的就是说明对称后的CGC同元CGC，是否整体差一个负号

        1，一些可以简化处理的情况可以走快速通道

        2，必须要遍历所有m的CGC的其他情况
        按照m_ν循环：
        2.1，整理6个矩阵中m_σ=1、m_μ=1、m_σ=h_σ、m_μ=h_μ的4个字典
        2.2，m_ν=1和m_ν=h_ν可以直接使用
        结束m_ν循环
        2.3，单σ=ν或μ=ν还是有简化算法可以榨取
        2.4，计算剩余4组ϵ

        6个矩阵指的是参数σ、μ、ν、τ下所有CGC，以m_σ、m_μ、m_ν为轴，构成一个三维tensor的6个表面
        它们分别代表固定m_σ=1、m_μ=1、m_ν=1、m_σ=h_σ、m_μ=h_μ、m_ν=h_ν下观察另外两个剩余参数

        注意：它所计算出的是书中所描述的ϵ。在多数情况下，ϵ(σμν)为正，它是适用的。但，当ϵ(σμν)为负时，就不适用了。
        所以，在实际使用中（指计算ISF和CGC的时候），我们给所有的ϵ再乘上一个ϵ(σμν)去使用，就没有问题了！
        """
        # TODO 所有简化情况重新检查，还不对，就暂时用最笨的方法实现，以后有需要再优化
        ϵ_dict = {}
        ϵ_flags = {}
        # 公共参数
        h_σ = data_sn.get_yt_num(data_σ_μ.σ)
        h_μ = data_sn.get_yt_num(data_σ_μ.μ)
        h_ν = data_sn.get_yt_num(ν)
        h_tuple = (h_σ, h_μ, h_ν)
        yd_σμν = (data_σ_μ.σ, data_σ_μ.μ, ν)
        Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in yd_σμν]

        if _debug_condition(data_σ_μ, ν=ν):
            logger.warning("######## catch meta data_σ_μ={}, ν={}".format(data_σ_μ, ν))

        # # 1，一些可以简化处理的情况可以走快速通道
        # if data_σ_μ.σ == data_σ_μ.μ == ν:
        #     # 这里，除了可以由σμν_0，σμν_1对称出所有σμν组合外，也不需要遍历所有m的CGC去凑一些dict了
        #     # 所以，先对m_ν=1、m_ν=h_ν做ϵ，再取对称就可以了
        #     # 1.1，m_ν=1
        #     flag, cgc_square_dict_1 = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, 1, is_flag_true_if_not_s_n=False)
        #     if not flag:
        #         err_msg = "get cgc_square_dict_1 with s_n={}, σ={}, μ={}, ν={}, τ={}, m_ν={} meet error with " \
        #                   "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, 1, cgc_square_dict_1)
        #         logger.error(err_msg)
        #         return False, err_msg
        #     adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict_1, data_σ_μ, ν, data_sn,
        #                                                        σμν_group=self.σμν_0, m_flags=(None, None, 1),
        #                                                        Λ_list_list=Λ_list_list)
        #     ϵ_dict.update(adding_ϵ_dict)
        #     ϵ_flags.update(adding_ϵ_flags)
        #
        #     # 1.2，m_ν=h_ν
        #     flag, cgc_square_dict_h = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, h_ν, is_flag_true_if_not_s_n=False)
        #     if not flag:
        #         err_msg = "get cgc_square_dict_h with s_n={}, σ={}, μ={}, ν={}, τ={}, h_ν={} meet error with " \
        #                   "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, h_ν, cgc_square_dict_h)
        #         logger.error(err_msg)
        #         return False, err_msg
        #     adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict_h, data_σ_μ, ν, data_sn,
        #                                                        σμν_group=self.σμν_1, m_flags=(None, None, h_ν),
        #                                                        Λ_list_list=Λ_list_list)
        #     ϵ_dict.update(adding_ϵ_dict)
        #     ϵ_flags.update(adding_ϵ_flags)
        #
        #     # 1.3，其他
        #     # 1.3.1 计算sym_ϵ_flags
        #     sym_ϵ_flags = {}
        #     for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
        #         # 对σμν_0和σμν_1采用两种对称（(1,3), (2,3)），实现全部ϵ
        #         s_13_d3, s_13_k4 = (d3[2], d3[1], d3[0]), (k4[2], k4[1], k4[0])
        #         s_23_d3, s_23_k4 = (d3[0], d3[2], d3[1]), (k4[0], k4[2], k4[1])
        #         s_13_ϵ_key = data_sn.get_ϵ_key_by_d3_k4(s_13_d3, s_13_k4)
        #         s_23_ϵ_key = data_sn.get_ϵ_key_by_d3_k4(s_23_d3, s_23_k4)
        #         sym_ϵ_flags[s_13_ϵ_key] = tuple(ϵ_flags[ϵ_key][d] if k is False else
        #                                         h_tuple[d] + 1 - ϵ_flags[ϵ_key][d] for d, k in zip(s_13_d3, s_13_k4))
        #         sym_ϵ_flags[s_23_ϵ_key] = tuple(ϵ_flags[ϵ_key][d] if k is False else
        #                                         h_tuple[d] + 1 - ϵ_flags[ϵ_key][d] for d, k in zip(s_23_d3, s_23_k4))
        #     # 1.3.2 整理一下sym_ϵ_flags
        #     m_to_load_dict = {}  # {m: [s_ϵ_key, ...]}
        #     for s_ϵ_key, s_flag in sym_ϵ_flags.items():
        #         if s_flag[2] not in m_to_load_dict:
        #             m_to_load_dict[s_flag[2]] = [s_ϵ_key]
        #         else:
        #             m_to_load_dict[s_flag[2]].append(s_ϵ_key)
        #
        #     # 1.3.3 根据整理后的sym_ϵ_flags找cgc
        #     adding_ϵ_dict = {}
        #     for m_to_load in m_to_load_dict.keys():
        #         if m_to_load == 1:
        #             cgc_square_dict_m = cgc_square_dict_1
        #         elif m_to_load == h_ν:
        #             cgc_square_dict_m = cgc_square_dict_h
        #         else:
        #             flag, cgc_square_dict_m = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_to_load,
        #                                                is_flag_true_if_not_s_n=False)
        #             if not flag:
        #                 err_msg = "get cgc_square_dict_m with s_n={}, σ={}, μ={}, ν={}, τ={}, m_to_load={} " \
        #                           "meet error with msg={}" \
        #                           "".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, h_ν, cgc_square_dict_m)
        #                 logger.error(err_msg)
        #                 return False, err_msg
        #         ϵ_key_list = m_to_load_dict[m_to_load]
        #         for ϵ_key in ϵ_key_list:
        #             ϵ_flag = sym_ϵ_flags[ϵ_key]
        #             cgc_key = (ϵ_flag[0], ϵ_flag[1])
        #             single_cgc_square = cgc_square_dict_m[cgc_key]  # 一定能拿到
        #             d3, k4 = data_sn.get_d3_k4_by_ϵ_key(ϵ_key)
        #             ΛΛ = 1
        #             for d, k in zip(d3, k4):
        #                 if k is True:
        #                     ΛΛ *= Λ_list_list[d][ϵ_flag[d] - 1]
        #             adding_ϵ_dict[ϵ_key] = sp.sign(single_cgc_square) * ΛΛ
        #     ϵ_dict.update(adding_ϵ_dict)
        #     ϵ_flags.update(sym_ϵ_flags)
        #
        #     # 1.4，简化情况可以提前返回
        #     if len(ϵ_dict) == len(ϵ_flags) == 24:
        #         return ϵ_dict, ϵ_flags
        #     else:
        #         err_msg = "len(ϵ_dict) == len(ϵ_flags) must eq 24 but not " \
        #                   "with ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags)
        #         return False, err_msg

        # 2，一般情况，必须要遍历CGC所有的m
        # 按照m_ν循环：
        # m_σ=1、m_μ=1、m_σ=h_σ、m_μ=h_μ的4个字典(不需要归一化，因为只需要它们的符号)
        m_μ_1_dict = {}
        m_μ_h_dict = {}
        m_σ_1_dict = {}
        m_σ_h_dict = {}
        for m_ν in range(1, h_ν + 1):
            flag, cgc_square_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_dict with s_n={}, σ={}, μ={}, ν={}, τ={}, m_ν={} meet error with " \
                          "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg

            if _debug_condition(data_σ_μ, ν=ν):
                logger.warning("## m_ν={}, cgc_square_dict={}".format(m_ν, cgc_square_dict))

            # 2.1，整理6个矩阵中m_σ=1、m_μ=1、m_σ=h_σ、m_μ=h_μ的4个字典
            for (m_σ, m_μ), v in cgc_square_dict.items():
                if m_μ == 1:
                    m_μ_1_dict[(m_σ, m_μ, m_ν)] = v
                if m_μ == h_μ:
                    m_μ_h_dict[(m_σ, m_μ, m_ν)] = v
                if m_σ == 1:
                    m_σ_1_dict[(m_σ, m_μ, m_ν)] = v
                if m_σ == h_σ:
                    m_σ_h_dict[(m_σ, m_μ, m_ν)] = v
            # 2.2，m_ν=1和m_ν=h_ν可以直接使用
            if m_ν == 1:
                adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict, data_σ_μ, ν, data_sn,
                                                                   σμν_group=self.σμν_0, m_flags=(None, None, m_ν),
                                                                   Λ_list_list=Λ_list_list)
                ϵ_dict.update(adding_ϵ_dict)
                ϵ_flags.update(adding_ϵ_flags)
            if m_ν == h_ν:
                adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict, data_σ_μ, ν, data_sn,
                                                                   σμν_group=self.σμν_1, m_flags=(None, None, m_ν),
                                                                   Λ_list_list=Λ_list_list)
                ϵ_dict.update(adding_ϵ_dict)
                ϵ_flags.update(adding_ϵ_flags)

            # if _debug_condition(data_σ_μ, ν=ν):
            #     logger.warning("#0# ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags))
        #
        # if _debug_condition(data_σ_μ, ν=ν):
        #     logger.warning("## m_μ_1_dict={}, m_μ_h_dict={}, m_σ_1_dict={}, m_σ_h_dict={}, "
        #                    "ϵ_dict={}, ϵ_flags={}"
        #                    "".format(m_μ_1_dict, m_μ_h_dict, m_σ_1_dict, m_σ_h_dict, ϵ_dict, ϵ_flags))

        # # TODO 这个简化有错误！它没考虑d的所有构型，先注释掉，如果性能提升大，以后再fix
        # # 2.3，单σ=ν或μ=ν还是有简化算法可以榨取
        # if data_σ_μ.σ == ν:
        #     for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
        #         # 对0和1采用σν对称可以算σμν_4、σμν_5
        #         if d3 == (0, 1, 2):  # 表示σμν（不看～）
        #             sym_d3, sym_k4 = (2, 1, 0), (k4[2], k4[1], k4[0])
        #         else:  # 这里必然是(1, 0, 2)  # 表示μσν（不看～）
        #             sym_d3, sym_k4 = (1, 2, 0), (k4[0], k4[2], k4[1])
        #         sym_ϵ_key = data_sn.get_ϵ_key_by_d3_k4(sym_d3, sym_k4)
        #         ϵ_dict[sym_ϵ_key], ϵ_flags[sym_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]
        #
        # if _debug_condition(data_σ_μ, ν=ν):
        #     logger.warning("#1# ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags))
        #
        # if data_σ_μ.μ == ν:
        #     for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
        #         # 对0和1采用μν对称可以算σμν_2、σμν_3
        #         if d3 == (0, 1, 2):  # 表示σμν（不看～）
        #             sym_d3, sym_k4 = (0, 2, 1), (k4[0], k4[2], k4[1])
        #         else:  # 这里必然是(1, 0, 2)  # 表示μσν（不看～）
        #             sym_d3, sym_k4 = (2, 0, 1), (k4[2], k4[1], k4[0])
        #         sym_ϵ_key = data_sn.get_ϵ_key_by_d3_k4(sym_d3, sym_k4)
        #         ϵ_dict[sym_ϵ_key], ϵ_flags[sym_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]
        #         if _debug_condition(data_σ_μ, ν=ν):
        #             logger.warning("## \nϵ_key={}, d3={}, k4={}, ϵ_dict[ϵ_key]={}, ϵ_flags[ϵ_key]={}\n"
        #                            "sym_ϵ_key={}, sym_d3={}, sym_k4={}, ϵ_dict[sym_ϵ_key]={}, ϵ_flags[sym_ϵ_key]={}"
        #                            "".format(ϵ_key, d3, k4, ϵ_dict[ϵ_key], ϵ_flags[ϵ_key],
        #                                      sym_ϵ_key, sym_d3, sym_k4, ϵ_dict[sym_ϵ_key], ϵ_flags[sym_ϵ_key]))

        # if _debug_condition(data_σ_μ, ν=ν):
        #     logger.warning("#2# ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags))

        # 2.4，计算剩余4组ϵ
        if any(i[0] not in ϵ_dict for i in self.σμν_2):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_μ_1_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_2, m_flags=(None, 1, None),
                                                               Λ_list_list=Λ_list_list)
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)
        if any(i[0] not in ϵ_dict for i in self.σμν_3):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_μ_h_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_3, m_flags=(None, h_μ, None),
                                                               Λ_list_list=Λ_list_list)
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)
        if any(i[0] not in ϵ_dict for i in self.σμν_4):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_σ_1_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_4, m_flags=(1, None, None),
                                                               Λ_list_list=Λ_list_list)
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)

            if _debug_condition(data_σ_μ, ν=ν):
                logger.warning("adding_ϵ_dict={}, adding_ϵ_flags={}".format(adding_ϵ_dict, adding_ϵ_flags))

        if any(i[0] not in ϵ_dict for i in self.σμν_5):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_σ_h_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_5, m_flags=(h_σ, None, None),
                                                               Λ_list_list=Λ_list_list)
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)

        # if _debug_condition(data_σ_μ, ν=ν):
        #     logger.warning("#3# ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags))

        if len(ϵ_dict) == len(ϵ_flags) == 24:
            return ϵ_dict, ϵ_flags
        else:
            err_msg = "len(ϵ_dict) == len(ϵ_flags) must eq 24 but not " \
                      "with ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags)
            return False, err_msg

    @staticmethod
    def calc_adding_ϵ(square_dict, data_σ_μ, ν, data_sn, σμν_group=None, m_flags=None, Λ_list_list=None):
        """根据传入的CGC，按照定义计算ϵ

        注意：中文版《群表示论的新途径》这段是错的，应该使用英文版《Group Representation Theory for Physicists》中的对应公式

        μσν:   对应ϵ1：公式4-116b
        νμσ:   对应ϵ2：公式4-118
        σνμ:   对应ϵ3：公式4-119
        σ~μ~ν: 对应ϵ4：公式4-122a
        σ~μν~: 对应ϵ5：公式4-122b
        σμ~ν~: 对应ϵ6：公式4-122c
        其他为自行推导，共计24个ϵ
        注意：square_dict不一定就是CGC字典，也可以是拼装出来的字典。如果是真CGC，它的key形如(σ,μ);如果是拼装的，key则形如(σ,μ,ν)
        σμν_group是预定的6个分组之一
        m_flags是len3tuple，非None值表示首先确定的杨盘的编号
        """
        adding_ϵ_dict = {}
        adding_ϵ_flags = {}

        # 公共条件
        default_σμν = (0, 1, 2)
        yd_σμν = (data_σ_μ.σ, data_σ_μ.μ, ν)
        if Λ_list_list is None:
            Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in yd_σμν]
        m_set_list = [set(j[i] for j in square_dict.keys()) if m_flags[i] is None else (m_flags[i],) for i in range(3)]

        # if _debug_condition(data_σ_μ, ν=ν):
        #     logger.warning("yd_σμν={}, Λ_list_list={}, m_set_list={}, square_dict={}, σμν_group={}, m_flags={}"
        #                    "".format(yd_σμν, Λ_list_list, m_set_list, square_dict, σμν_group, m_flags))

        for ϵ_key, d3, k4 in σμν_group:
            # ϵ_key是字符串形式的键，p.s. μσ~ν~； d3是组合的数字表示，p.s. (1, 0, 2)； k4是组合的共轭bool，p.s.

            # if _debug_condition(data_σ_μ, ν=ν):
            #     logger.warning("ϵ_key={}, d3={}, k4={}, adding_ϵ_dict={}"
            #                    "".format(ϵ_key, d3, k4, adding_ϵ_dict))

            if ϵ_key in adding_ϵ_dict:
                continue

            # if _debug_condition(data_σ_μ, ν=ν) and ϵ_key == "μν~σ~":
            #     logger.warning("ϵ_key={}, d3={}, k4={}".format(ϵ_key, d3, k4))

            # 下面按顺位缩小范围
            # 第一顺位
            first_σμν_index = d3[0]
            first_σμν_set = m_set_list[first_σμν_index]
            first_m = min(first_σμν_set) if k4[0] is False else max(first_σμν_set)

            # 第二顺位
            second_σμν_index = d3[1]
            # 无论是len2key还是len3key，它们都可以使用下面的处理
            second_σμν_set = set(i[second_σμν_index] for i in square_dict.keys() if i[first_σμν_index] == first_m)
            second_m = min(second_σμν_set) if k4[1] is False else max(second_σμν_set)

            # 因为第三顺位是最先确定的条件，所以现在可以拿到m_tuple了，它表示例如σμν的非首项坐标m_tuple在做对称后将成为μσν的CGC的首项
            m_tmp_list = list(m_flags)
            m_tmp_list[first_σμν_index] = first_m
            m_tmp_list[second_σμν_index] = second_m
            m_tuple = tuple(m_tmp_list)  # 是meta的坐标哦！！！

            # if _debug_condition(data_σ_μ, ν=ν) and ϵ_key == "μν~σ~":
            #     logger.warning("ϵ_key={}, m_tuple={}".format(ϵ_key, m_tuple))

            # 根据前面计算Λ * Λ
            ΛΛ = 1
            for d, k in zip(d3, k4):
                if k is True:
                    ΛΛ *= Λ_list_list[d][m_tuple[d] - 1]

            square_dict_key = m_tuple[:-1:] if d3[-1] == default_σμν[-1] else m_tuple  # 真CGC是len2key，拼装的len3key
            adding_ϵ_dict[ϵ_key] = sp.sign(square_dict[square_dict_key]) * ΛΛ
            adding_ϵ_flags[ϵ_key] = m_tuple

            # if _debug_condition(data_σ_μ, ν=ν) and ϵ_key == "μν~σ~":
            #     logger.warning("adding_ϵ_dict[ϵ_key]={}, square_dict_key={}, square_dict[square_dict_key]={}, ΛΛ={}"
            #                    "".format(adding_ϵ_dict[ϵ_key], square_dict_key, square_dict[square_dict_key], ΛΛ))

            # # 一二顺位相等时的简化计算
            # # TODO 这个简化有错误！它没考虑d的所有构型，先注释掉，如果性能提升大，以后再fix
            # if yd_σμν[first_σμν_index] == yd_σμν[second_σμν_index]:
            #     sym_d3, sym_k4 = (d3[1], d3[0], d3[2]), (k4[1], k4[0], k4[2])
            #     sym_ϵ_key = data_sn.get_ϵ_key_by_d3_k4(sym_d3, sym_k4)
            #     adding_ϵ_dict[sym_ϵ_key] = adding_ϵ_dict[ϵ_key]
            #     adding_ϵ_flags[sym_ϵ_key] = adding_ϵ_flags[ϵ_key]
            #
            #     if _debug_condition(data_σ_μ, ν=ν):
            #         logger.warning("sym_ϵ_key={}, adding_ϵ_dict[sym_ϵ_key]={}, adding_ϵ_dict[ϵ_key]={}\n"
            #                        "first_σμν_index={}, second_σμν_index={}"
            #                        "".format(sym_ϵ_key, adding_ϵ_dict[sym_ϵ_key], adding_ϵ_dict[ϵ_key],
            #                                  first_σμν_index, second_σμν_index))

        return adding_ϵ_dict, adding_ϵ_flags

    def _calc_single_sym_mode_ϵ_by_meta(self, sym_mode, sym_mode_d3, sym_mode_k4, meta_ϵ_dict, meta_ϵ_flags,
                                        meta_Λ_list_list, sym_Λ_list_list, meta_h_list, data_sn):
        """
        根据指定的sym_mode求其24个sym_key

        注意，下面的代码中将出现五套四个σμν，分别表示：
        meta_mode：就是meta(σμν)！
        meta_key：形如ν~σμ~：表示meta(σμν)有24个ϵ，每个ϵ的具体对称形式和键就是meta_key
        sym_mode：形如ν~μσ~(=σ'μ'ν')：表示meta(σμν)有24种对称，当前把meta变换成sym_mode，以sym_mode为原点，sym_mode也将有24个自己的ϵ
        sym_mode'：形如σ'μ'ν'(=ν~μσ~)：以sym_mode作为原点的记法
        sym_key：形如σ'ν'~μ'~：表示sym_mode也将有24个自己的ϵ，每个ϵ的具体对称形式和键就是sym_key
        上面举例中，meta的ν~σμ~和sym的σ'ν'~μ'~对应的是同一个实体

        ΛΛ^xxx(m|yyy)：表示按照xxx取Λ函数，按照m取函数的值。代码中出现三种：
        形如ΛΛ_d^σ'ν'~μ'~(m'| sym_key cond)：表示定义式Λ^μ'(h_μ') * Λ^ν'(max(m_ν')|min(m_σ'))
        形如ΛΛ_u^ν~μσ~(m｜ν~σμ~ cond)：表示使用式Λ^ν(max(m_ν)) * Λ^σ(min(m_σ))

        请注意定义式与使用式关于角标的不同，特别是ΛΛ部分
        定义式：ϵ(meta_key) = sign(CGC^meta_mode[m｜meta_key cond]) * ΛΛ_d^meta_key(m｜meta_key cond)
        使用式：CGC^sym_mode(m|meta_key cond)
            = ϵ(σμν) * h_σ/h_ν * ϵ(sym_mode) * ΛΛ_u^sym_mode(m｜meta_key cond) * CGC^meta_mode[m｜meta_key cond]

        从结果上看，这里所求的是一种确定的对称组合sym_mode'(σ'μ'ν')下，它的24个ϵ
        其单独的一个ϵ的组合是sym_key(σ'ν'~μ'~)，那么按照定义，其ϵ可表示为两部分
        sign(sym_mode'经过sym_key对称后，是sym_key原点的坐标的sym_mode'的CGC) * ΛΛ_d^sym_key(sym_key相对于sym_mode的～关系)
        记为：sign(CGC^σ'μ'ν'[m'| sym_key cond]) * ΛΛ_d^σ'ν'~μ'~(m'| sym_key cond)
        首先是m'| sym_key cond
        此m'是sym_key相对于sym_mode'的
        因为sym_mode也是可以被meta_mode+meta_key表示的
        我们也可以根据meta_key将它与meta_mode对应起来，既：m｜meta_key(m'| sym_key cond)
        注意，这里不需要知道m'的具体值，只需要知道其规则，就可以反推了

        例如：
        meta_mode(σμν) = [4, 1] [3, 2] [3, 1, 1]
        存在ϵ(ν~μσ~)
        则有
        sym_mode'(σ'μ'ν') = [3, 1, 1] [3, 2] [2, 1, 1, 1] (对于meta，它可以通过ϵ(ν~μσ~)得到)
        现在求sym_key=ϵ'(σ'ν'~μ'~)的m'（既[3, 1, 1], [4, 1], [2, 2, 1]）
        按照定义，h_μ'=h(m_[3, 2]')=5；min(m_σ')=min(m_[3, 1, 1]')=;max(m_ν')=max(m_[2, 1, 1, 1]')=
        上述m'不能通过未知的当前CGC确定！
        但是，可以通过sym_mode(ϵ(ν~μσ~))将m'对σ'μ'ν'的关系转化为m对σμν的关系
        既，(μ'=μ)(h_μ')=h_μ=h(m_[3, 2])=5;(σ'=ν~)(min(m_σ'))=max(m_ν)=;(ν'=σ~)(max(m_ν'))=min(m_σ)=
        这里，同样也不需要真正去查表得到m，只需要看24个meta_key谁与之对应即可
        得到meta_key(ν~σμ~)满足条件，
        并且，它是完全已知的(m_μ=5;m_ν=6;m_σ=3;CGC^σμν(m|ν~σμ~)=-3/8;h_ν=6;h_μ=5;h_σ=4)  # ϵ(ν~σμ~)=sign(-3/8*-*+)=1;
        为了验证，我们正向再走一遍：
        将上述m(m_σ=3;m_μ=5;m_ν=6)按照sym_mode=ϵ(ν~μσ~)变换(m_ν~=1;m_μ=5;m_σ~=2)，得到CGC^ν~μσ~，它也是CGC^σ'μ'ν'，
        有：
        CGC^ν~μσ~(m|ν~σμ~)
        # 注意，这里不是 定义式 而是 使用式，所以h/h和ϵ用sym_mode 而CGC用meta_key condition
        # 但是ΛΛ最特殊，在σμν的选取上，它使用sym_mode；而在m的取值上，它使用meta_key !!!!
        = ϵ(σμν) * h_σ/h_ν * ϵ(sym_mode) * ΛΛ_u^sym_mode(m｜meta_key cond) * CGC^meta_mode[m｜meta_key cond]
        = ϵ(σμν) * h_σ/h_ν * ϵ(ν~μσ~) * ΛΛ_u^ν~μσ~(m｜ν~σμ~ cond) * CGC^σμν[m｜ν~σμ~ cond]
        = ϵ(σμν) * h_σ/h_ν * ϵ(ν~μσ~) * ΛσΛν(m_σ=3;m_ν=6) * CGC^σμν[m_σ=3;m_μ=5;m_ν=6]
        # ϵ(ν~μσ~)=ϵ(m_σ=4;m_ν=6;m_μ=2)=sign(CGC^σμν(m|ν~μσ~)*Λ_σ(4)*Λ_ν(6))=sign(1/8*-*+)=-1
        # Λ_σ(3)=+;Λ_ν(6)=+
        = + * 4/6 * -1 * +*+ * -3/8 = 1/4
        它就是正向使用meta_mode+meta_key得到sym_mode'的对称公式结果，查表验证确实等于：
        因为对应关系m_σ'=(h_ν+1-m_ν)=6+1-6=1; m_μ'=m_μ=5; m_ν'=(h_σ+1-m_σ)=4+1-3=2
        CGC^σ'μ'ν'[m_σ'=1; m_μ'=5; m_ν'=2] = 1/4
        以(σ'μ'ν')为原点，把当前m'，经过sym_key(σ'=[3, 1, 1];ν'~=[4, 1];μ'~=[2, 2, 1])对称后，变为(m_σ'=1;m_ν'~=3;m_μ'~=1)
        查表可知，它就是σ'=[3, 1, 1];ν'~=[4, 1];μ'~=[2, 2, 1]的首项。正向是正确的！
        通过上述方式，我们可以确定m'| sym_key condition，同时，还找到了这个CGC在meta_mode中的对应
        到此为止，ΛΛ_d^sym_key(m'| sym_key cond)已经完全确定。

        下面仍以此例子，展示sign(CGC^sym_mode'[m'| sym_key cond])的推导过程：
        首先，根据m'副产品可知，
        CGC^σ'μ'ν'[m'| sym_key cond]
        = ϵ(σμν) * h_σ/h_ν * ϵ(sym_mode) * ΛΛ_u^ν~μσ~(m｜ν~σμ~ cond) * CGC^σμν[m｜meta_key cond]
        注意：此处比书中公式多了一个ϵ(σμν)，是为了满足τ>1时，出现的首项为负的情况。如果首项为正，则可以省略
        套上sign后，可简化为：ϵ(σμν) * ϵ(sym_mode) * ΛΛ_u^sym_mode(m｜meta_key cond) * sign(CGC^σμν[m｜meta_key cond])
        这里，前三项都是易知的，问题就简化为求sign(CGC^σμν[m｜meta_key cond])
        在例子中，就是求sign(CGC^σμν[m_σ=3;m_μ=5;m_ν=6])
        我们这里不希望使用真正的CGC，因为只需求sign，使用ϵ的定义式就足够了
        由ϵ的定义式：ϵ(meta_key) = sign(CGC^meta_mode[m｜meta_key cond]) * ΛΛ_d^meta_key(m｜meta_key cond)
        有：ϵ(ν~σμ~) = sign(CGC^σμν[m｜ν~σμ~]) * ΛΛ_d^ν~σμ~(m｜ν~σμ~)
        所以，sign(CGC^σμν[m｜ν~σμ~]) = ϵ(ν~σμ~) / ΛΛ_d^ν~σμ~(m｜ν~σμ~) = ϵ(ν~σμ~) * ΛΛ_d^ν~σμ~(m｜ν~σμ~)
        例如，
        sign(CGC^σμν[m_σ=3;m_μ=5;m_ν=6])
        = ϵ(ν~σμ~) * ΛΛ_d^ν~σμ~(m_σ=3;m_μ=5;m_ν=6)
        = ϵ(ν~σμ~) * Λ_ν(m_ν=6) * Λ_μ(m_μ=5)

        至此，sym_mode的ϵ(sym_key)，已经可以由已知meta_mode+计算出的对应ϵ(meta_key)表示，记为：
        ϵ(sym_key)
        = sign(CGC^sym_mode'[m'| sym_key cond]) * ΛΛ_d^sym_key(m'| sym_key cond)
        = ϵ(meta_mode) * ϵ(sym_mode) * sign(CGC^meta_mode[m｜meta_key cond]) * ΛΛ_u^sym_mode(m｜meta_key cond)
          * ΛΛ_d^sym_key(m'| sym_key cond)
        = ϵ(meta_mode) * ϵ(sym_mode)
          * ϵ(meta_key) * ΛΛ_d^meta_key(m｜meta_key cond)
          * ΛΛ_u^sym_mode(m｜meta_key cond) * ΛΛ_d^sym_key(m'| sym_key cond)
        = ϵ(meta_mode) * ϵ(sym_mode) * ϵ(meta_key) * ΛΛ_d^meta_key(m｜meta_key cond)
                                                   * ΛΛ_u^sym_mode(m｜meta_key cond)
                                                   * ΛΛ_d^sym_key(m'| sym_key cond)
        例如：
        ϵ(σ'ν'~μ'~)  # sym_mode=ν~μσ~
        首先计算出meta_key=ν~σμ~
        = ϵ(σμν) * ϵ(ν~μσ~) * ϵ(ν~σμ~) * ΛΛ_d^ν~σμ~(m_σ=3;m_μ=5;m_ν=6)
                                       * ΛΛ_u^ν~μσ~(m_σ=3;m_μ=5;m_ν=6; m_σ'=1;m_μ'=5;m_ν'=2)
                                       * ΛΛ_d^σ'ν'~μ'~(m_σ'=1;m_μ'=5;m_ν'=2)
        = ϵ(σμν) * ϵ(ν~μσ~) * ϵ(ν~σμ~) * Λ_ν(m_ν=6) * Λ_μ(m_μ=5)
                                       * Λ_ν(m_ν=6) * Λ_σ(m_σ=3)
                                       * Λ_ν'(m_ν'=2) * Λ_μ'(m_μ'=5)
        """
        sym_ϵ_dict = {}
        sym_ϵ_flags = {}
        default_meta_mode = "σμν"
        # 1, ϵ(meta_mode)
        ϵ_meta_mode = meta_ϵ_dict.get(default_meta_mode)
        # 2, ϵ(sym_mode)
        ϵ_sym_mode = meta_ϵ_dict.get(sym_mode)
        for sym_key, sym_key_d3, sym_key_k4 in \
                chain(self.σμν_0, self.σμν_1, self.σμν_2, self.σμν_3, self.σμν_4, self.σμν_5):
            # 求meta_key
            meta_key_d3 = tuple(sym_mode_d3[i] for i in sym_key_d3)
            meta_key_k4 = tuple(sym_mode_k4[i] is not k for i, k in zip(sym_key_d3, sym_key_k4))
            meta_key = data_sn.get_ϵ_key_by_d3_k4(meta_key_d3, meta_key_k4)
            # 3, ϵ(meta_key)
            ϵ_meta_key = meta_ϵ_dict.get(meta_key)
            # 4, ΛΛ_d^meta_key(m｜meta_key cond)
            m_meta_key_cond = meta_ϵ_flags.get(meta_key)  # σμν顺序
            ΛΛ_d_meta_key = 1
            for d, k in zip(meta_key_d3, meta_key_k4):
                if k is True:
                    ΛΛ_d_meta_key *= meta_Λ_list_list[d][m_meta_key_cond[d] - 1]
            # 5, ΛΛ_u^sym_mode(m｜meta_key cond)
            ΛΛ_u_sym_mode = 1
            for d, k in zip(sym_mode_d3, sym_mode_k4):
                if k is True:
                    ΛΛ_u_sym_mode *= meta_Λ_list_list[d][m_meta_key_cond[d] - 1]
            # 6, ΛΛ_d^sym_key(m'| sym_key cond)
            m_sym_key_cond = tuple(m_meta_key_cond[i] if j is False else meta_h_list[i] + 1 - m_meta_key_cond[i]
                                   for i, j in zip(sym_mode_d3, sym_mode_k4))  # σ'μ'ν'顺序
            ΛΛ_d_sym_key = 1
            for d, k in zip(sym_key_d3, sym_key_k4):
                if k is True:
                    ΛΛ_d_sym_key *= sym_Λ_list_list[d][m_sym_key_cond[d] - 1]
            # 总结
            sym_ϵ_dict[sym_key] = ϵ_meta_mode * ϵ_sym_mode * ϵ_meta_key * ΛΛ_d_meta_key * ΛΛ_u_sym_mode * ΛΛ_d_sym_key
            sym_ϵ_flags[sym_key] = m_sym_key_cond

        return sym_ϵ_dict, sym_ϵ_flags

    def calc_symmetry_ϵ_by_meta_include_save(self, meta_ϵ_dict, meta_ϵ_flags, meta_data_σ_μ, meta_ν, τ, data_sn):
        """
        根据元ϵ计算对称的ϵ（坐标系的原点变换）
        """
        meta_yd_σμν = (meta_data_σ_μ.σ, meta_data_σ_μ.μ, meta_ν)
        meta_Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in meta_yd_σμν]
        meta_h_list = [data_sn.get_yt_num(yd) for yd in meta_yd_σμν]
        for sym_mode, sym_mode_d3, sym_mode_k4 in \
                chain(self.σμν_0, self.σμν_1, self.σμν_2, self.σμν_3, self.σμν_4, self.σμν_5):

            # if _debug_condition(meta_data_σ_μ):
            #     logger.warning("@@@@ meta_data_σ_μ={}, "
            #                    "sym_mode={}, sym_mode_d3={}, sym_mode_k4={}, "
            #                    "data_sn.yd_tilde_list={}".format(meta_data_σ_μ, sym_mode, sym_mode_d3, sym_mode_k4,
            #                                                      data_sn.yd_tilde_list))

            sym_yd_σμν = tuple(meta_yd_σμν[d] if k is False else
                               data_sn.get_tilde(meta_yd_σμν[d]) for d, k in zip(sym_mode_d3, sym_mode_k4))

            # if _debug_condition(meta_data_σ_μ):
            #     logger.warning("@@@@ meta_data_σ_μ={}".format(meta_data_σ_μ))

            sym_ϵ_tuple = (self.s_n, *sym_yd_σμν, τ)
            _, is_calc_ed = is_ϵ_exist(*sym_ϵ_tuple)
            if is_calc_ed is True:
                continue
            sym_Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in sym_yd_σμν]

            # if _debug_condition(meta_data_σ_μ):
            #     logger.warning("@@@@ meta_data_σ_μ={}, meta_yd_σμν={}, sym_yd_σμν={}, "
            #                    "meta_Λ_list_list={}, sym_Λ_list_list={}".format(meta_data_σ_μ, meta_yd_σμν,
            #                                                                     tuple(sym_yd_σμν),
            #                                                                     meta_Λ_list_list, sym_Λ_list_list))

            sym_ϵ_dict, sym_ϵ_flags = \
                self._calc_single_sym_mode_ϵ_by_meta(sym_mode, sym_mode_d3, sym_mode_k4, meta_ϵ_dict, meta_ϵ_flags,
                                                     meta_Λ_list_list, sym_Λ_list_list, meta_h_list, data_sn)
            if len(sym_ϵ_dict) == len(sym_ϵ_flags) == 24:
                flag, msg = save_ϵ(*(*sym_ϵ_tuple, sym_ϵ_dict, sym_ϵ_flags))
                if not flag:
                    err_msg = "save_ϵ fail by sym_ϵ_tuple={}, sym_ϵ_dict={}, sym_ϵ_flags={} with " \
                              "msg={}".format(sym_ϵ_tuple, sym_ϵ_dict, sym_ϵ_flags, msg)
                    logger.error(err_msg)
                    return False, err_msg
            else:
                err_msg = "_calc_single_sym_mode_ϵ_by_meta return wrong sym_ϵ_dict={}, sym_ϵ_flags={}".format(
                    sym_ϵ_tuple, sym_ϵ_dict, sym_ϵ_flags)
                logger.error(err_msg)
                return False, err_msg

        return True, None
