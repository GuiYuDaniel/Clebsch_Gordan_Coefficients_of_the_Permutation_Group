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
from itertools import product, combinations, combinations_with_replacement, chain
from conf.cgc_config import default_s_n, group_d3, group_k4
from conf.cgc_config import min_s_n_of_isf, min_s_n_of_cgc, min_s_n_of_branching_law, min_s_n_of_ϵ
from core.young_diagrams import load_young_diagrams, calc_young_diagram_tilde
from core.young_tableaux import load_young_table_num, load_young_table_phase_factor
from core.branching_laws import load_branching_law
from core.yamanouchi_matrix import load_yamanouchi_matrix
from core.cg_series import load_cg_series
from core.eigenvalues import load_eigenvalues
from core.cgc_utils.cgc_db_typing import ISFInfo, CGCInfo, EInfo
from core.cgc_utils.cgc_local_db import get_isf_file_name, get_isf_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_cgc_file_name, get_cgc_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_ϵ_file_name, get_ϵ_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


def _debug_condition(data_σ_μ, ν_st):
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
    if cond_21 or cond_23:
        return True
    return False


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
    isf_func = ISFHelper()
    cgc_func = CGCHelper()
    ϵ_func = EHelper()
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
        data_st = data_si if data_si is not None else DataHelper(s_i - 1)
        # 对于Sn=1、2，是不启动add_more_data的，后面的函数也不会调用它
        if data_si is None and s_i - 1 - 1 >= min_s_n_of_branching_law:
            data_st_st_yt_num_list = []
            _, yd_st_st_list = load_young_diagrams(s_i - 1 - 1)
            for yd_st_st in yd_st_st_list:
                _, total_num = load_young_table_num(s_i - 1 - 1, yd_st_st)
                data_st_st_yt_num_list.append(total_num)
            data_st.add_more_data(data_st_st_yt_num_list, yd_st_st_list)
        # 对于Sn，需要计算本轮的准备数据
        data_si = DataHelper(s_i)
        data_si.add_more_data(data_st.yt_num_list, data_st.yd_list)

        # σ μ
        '''
        第六版算法如下：

        生成一组σμ组合：
        
        1，判断，跳过非元σμ组合，并初始化辅助数据
        
        2，计算元σμ下所有ISF、CGC和ϵ（使用4-195a式而不是ϵ）
        ν'循环：
        2.1，计算row_index_tmp_list，并判断ν'的合法性
        2.2，计算σμν'的ISF
        2.2.1，使用4-189a和4-195a式计算ν'对应的ISF¥[σ][μ][ν']¥
        2.2.2，保存和注册元ISF信息和数据
        2.3，根据当前ISF计算元CGC以及ϵ【[σ][μ][ν]τ】
        2.3.1，计算元CGC
        2.3.2，计算8个ϵ
        注意：不用管什么对称性，直接按照定义计算上述所有参数即可
        
        3，计算σμ外对称ISF
        metaISF循环：(元ISF的标志是ν')
        3.0，数据准备
        3.1，根据ϵ1算¥[μ][σ][ν']¥
        3.2，根据ϵ4算¥[σ~][μ~][ν‘]¥
        3.3，根据ϵ14算¥[μ~][σ~][ν’]¥
        3.4，根据ϵ5算¥[σ~][μ][ν‘~]¥
        3.5，根据ϵ15算¥[μ][σ~][ν’~]¥
        3.6，根据ϵ6算¥[σ][μ~][ν‘~]¥
        3.7，根据ϵ16算¥[μ~][σ][ν’~]¥
        
        4，计算σμ外对称ϵ
        metaCGC循环：(元ϵ的标志是ν+τ)
        4.0，数据准备
        4.1，算【[μ][σ][ν]τ】ϵ1_dict
        4.2，算【[σ~][μ~][ν]τ】ϵ4_dict
        4.3，算【[μ~][σ~][ν]τ】ϵ14_dict
        4.4，算【[σ~][μ][ν~]τ】ϵ5_dict
        4.5，算【[μ][σ~][ν~]τ】ϵ15_dict
        4.6，算【[σ][μ~][ν~]τ】ϵ6_dict
        4.7，算【[μ~][σ][ν~]τ】ϵ16_dict
        
        5，计算σμ外对称CGC
        metaCGC循环：(元CGC的标志是ν+τ)
        5.0，数据准备
        5.1，根据ϵ1算【[μ][σ][ν]τ】
        5.2，根据ϵ4算【[σ~][μ~][ν]τ】
        5.3，根据ϵ14算【[μ~][σ~][ν]τ】
        5.4，根据ϵ5算【[σ~][μ][ν~]τ】
        5.5，根据ϵ15算【[μ][σ~][ν~]τ】
        5.6，根据ϵ6算【[σ][μ~][ν~]τ】
        5.7，根据ϵ16算【[μ~][σ][ν~]τ】
        '''
        for σ, μ in combinations_with_replacement(data_si.yd_list, 2):  # [σ], [μ]双循环，组合不排列
            # 1，判断，跳过非元σμ组合，并初始化辅助数据
            _, is_σ_μ_calc_ed = isf_func.is_isf_σ_μ_exist(σ, μ)
            if is_σ_μ_calc_ed is True:
                continue

            logger.debug("meta σ={}, μ={}".format(σ, μ))

            # 初始化元σμ循环数据
            data_σ_μ = ΣMDataHelper(s_i, σ, μ, data_si, data_st)

            # 2，计算元σμ下所有ISF、CGC和ϵ（使用4-195a式而不是ϵ）
            for ν_st in data_st.yd_list:  # ν_st比ν在前

                # 2.1，计算row_index_tmp_list，并判断ν'的合法性
                '''
                所有情况枚举：
                a，σμ=μσ，已经被combinations_with_replacement排除了
                b，is_σ_μ_calc_ed为真，说明本σμ，已经在元σμ组合中，被对称性计算过了
                以下情况，是本组合内，需要判断的。同时，本组合必是元σμ组合
                c，row_index_tmp_list=[]的，表示该组合下不存在ISF（跳过）
                d，其他
                '''
                # c，row_index_tmp_list=[]的，表示该组合下不存在ISF（跳过）
                row_index_tmp_list = isf_func.calc_row_indexes_tmp(ν_st, data_st, data_σ_μ)  # 带着None的rows
                if not row_index_tmp_list:
                    continue
                logger.debug("new meta combination σ={}, μ={}, ν_st={}".format(σ, μ, ν_st))

                # d & 2.2，计算σμν'的ISF
                '''
                isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                                   "cols":[[ν], ([ν], τ), ...],  # 有自由度tuple，无自由度list
                                   "isf": isf_square_matrix}  # sp.Matrix([len(rows), len(cols)])
                '''
                # 2.2.1，使用4-189a和4-195a式计算ν'对应的ISF¥[σ][μ][ν']¥
                meta_isf_start_time = time.time()
                flag, meta_isf_square_dict = isf_func.calc_meta_isf_dict(row_index_tmp_list, ν_st,
                                                                         data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_meta_isf_dict meet error by row_index_tmp_list={}, ν_st={}, data_si={}, " \
                              "data_st={}, data_σ_μ={} with msg={}".format(row_index_tmp_list, ν_st, data_si,
                                                                           data_st, data_σ_μ, meta_isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
                # 2.2.2，保存和注册元ISF信息和数据
                data_σ_μ.register_isf_ν_st_list(ν_st)
                data_σ_μ.register_isf_square_dict_by_ν_st(ν_st, meta_isf_square_dict)
                meta_isf_speed_time = int(time.time() - meta_isf_start_time)
                isf_speed_time_s_i += meta_isf_speed_time
                flag, msg = save_isf(s_i, σ, μ, ν_st, meta_isf_square_dict, meta_isf_speed_time)
                if not flag:
                    err_msg = "save meta isf meet error with s_i={}, σ={}, μ={}, ν_st={}, meta_isf_square_dict={}, " \
                              "msg={}".format(s_i, σ, μ, ν_st, meta_isf_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

                if _debug_condition(data_σ_μ, ν_st):
                    logger.warning("@@@@ meta_isf_square_dict={}".format(meta_isf_square_dict))

                # 2.3，根据当前ISF计算元CGC【[σ][μ][ν]τ】
                meta_cgc_start_time = time.time()
                flag, msg = cgc_func.calc_meta_cgc_include_save(meta_isf_square_dict, ν_st, data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_meta_cgc_include_save fail by meta_isf_square_dict={}, ν_st={}, " \
                              "data_si={}, data_st={}, data_σ_μ={} with " \
                              "msg={}".format(meta_isf_square_dict, ν_st, data_si, data_st, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg
                meta_cgc_speed_time = int(time.time() - meta_cgc_start_time)
                cgc_speed_time_s_i += meta_cgc_speed_time

            meta_ϵ_start_time = time.time()
            meta_ϵ_ν_τ_list = data_σ_μ.get_meta_cgc_ν_τ_list()  # 这里借用cgc的就行
            for (ν, τ) in meta_ϵ_ν_τ_list:
                # 3，计算当前σμ组合全部ν的24个元ϵ
                ϵ_dict, ϵ_flags = ϵ_func.calc_meta_ϵ_dict(data_σ_μ, ν, τ, data_si)
                if not ϵ_dict:
                    err_msg = "calc_meta_ϵ_dict fail by data_σ_μ={}, ν={}, τ={}, data_si={} " \
                              "with msg={}".format(data_σ_μ, ν, τ, data_si, ϵ_flags)
                    logger.error(err_msg)
                    return False, err_msg
                flag, msg = save_ϵ(s_i, data_σ_μ.σ, data_σ_μ.μ, ν, τ, ϵ_dict, ϵ_flags)
                if not flag:
                    err_msg = "save_ϵ fail by s_i={}, σ={}, μ={}, ν={}, τ={}, ϵ_dict={}, ϵ_flags={} with " \
                              "msg={}".format(s_i, data_σ_μ.σ, data_σ_μ.μ, ν, τ, ϵ_dict, ϵ_flags, msg)
                    logger.error(err_msg)
                    return False, err_msg
                # 4，计算剩余23组非元对称组合的自身的24个ϵ

                flag, msg = cgc_func.calc_symmetry_ϵ_by_meta_include_save(ν, τ, data_si, data_σ_μ)
                if not flag:
                    err_msg = "calc_symmetry_ϵ_by_meta_include_save fail by ν={}, τ={}, data_si={}, data_σ_μ={} " \
                              "with msg={}".format(ν, τ, data_si, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg

            meta_ϵ_speed_time = int(time.time() - meta_ϵ_start_time)
            ϵ_speed_time_s_i += meta_ϵ_speed_time

            # 4，计算剩余23组非元对称组合的自身的24个ϵ

            # 4，计算全部σμ外对称ISF
            symmetry_isf_start_time = time.time()
            meta_isf_ν_st_list = data_σ_μ.get_isf_ν_st_list()
            logger.debug("meta_isf_ν_st_list={}".format(meta_isf_ν_st_list))
            for meta_ν_st in meta_isf_ν_st_list:
                meta_isf_square_dict = data_σ_μ.get_isf_square_dict_by_ν_st(meta_ν_st)
                flag, msg = isf_func.calc_symmetry_isf_include_save(meta_isf_square_dict, meta_ν_st,
                                                                    data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_symmetry_isf_include_save fail by meta_isf_square_dict={}, meta_ν_st={}, " \
                              "data_sn={}, data_st={}, data_σ_μ={} with " \
                              "msg={}".format(meta_isf_square_dict, meta_ν_st, data_si, data_st, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg
            symmetry_isf_speed_time = int(time.time() - symmetry_isf_start_time)
            isf_speed_time_s_i += symmetry_isf_speed_time

            # 4，计算σμ外对称ϵ
            symmetry_ϵ_start_time = time.time()
            meta_ϵ_ν_τ_list = data_σ_μ.get_meta_cgc_ν_τ_list()  # 这里借用cgc的就行
            for (ν, τ) in meta_ϵ_ν_τ_list:
                flag, msg = cgc_func.calc_symmetry_ϵ_by_meta_include_save(ν, τ, data_si, data_σ_μ)
                if not flag:
                    err_msg = "calc_symmetry_ϵ_by_meta_include_save fail by ν={}, τ={}, data_si={}, data_σ_μ={} " \
                              "with msg={}".format(ν, τ, data_si, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg
            symmetry_ϵ_speed_time = int(time.time() - symmetry_ϵ_start_time)
            ϵ_speed_time_s_i += symmetry_ϵ_speed_time

            # 5，计算σμ外对称CGC
            symmetry_cgc_start_time = time.time()
            meta_cgc_ν_τ_list = data_σ_μ.get_meta_cgc_ν_τ_list()
            logger.debug("meta_cgc_ν_τ_list={}".format(meta_cgc_ν_τ_list))
            for (ν, τ) in meta_cgc_ν_τ_list:
                flag, msg = cgc_func.calc_symmetry_cgc_include_save(ν, τ, data_si, data_σ_μ)
                if not flag:
                    err_msg = "calc_symmetry_cgc_include_save fail by ν={}, τ={}, data_si={}, data_σ_μ={} with " \
                              "msg={}".format(ν, τ, data_si, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg
            symmetry_cgc_speed_time = int(time.time() - symmetry_cgc_start_time)
            cgc_speed_time_s_i += symmetry_cgc_speed_time

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


def save_isf(s_n: int, σ: list, μ: list, ν_st: list, isf_square_dict: dict, speed_time: int):
    """
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

    db_info = ISFInfo(s_n)
    _, file_name = get_isf_file_name(s_n, σ, μ, ν_st)
    table = {"file_name": file_name,
             "data": isf_square_dict,
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
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
    else:
        err_msg = "output_mode={} must in ['', 'single_row', 'single_col', 'single_isf'] but not".format(output_mode)
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
        isf_square_dict = data.get("data")
        if isinstance(isf_square_dict, dict) and isf_square_dict:
            if output_mode == "single_row":
                # 输出单行，额外输入行指标([σ'], [μ'], τ') or ([σ'], [μ'])
                ex_row = ex_params[: -1] if None in ex_params else ex_params  # 使它也能接受([σ'], [μ'], None)
                if ex_row in isf_square_dict.get("rows", []):
                    row_index = isf_square_dict.get("rows").index(ex_row)
                    rst_dict = {"cols": isf_square_dict["cols"],
                                "single_row": isf_square_dict["isf"][row_index, :]}
                    return True, rst_dict  # bingo(1/4)！
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
                    return True, rst_dict  # bingo(2/4)！
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
                    return True, rst  # bingo(3/4)！
                else:
                    err_msg = "ex_row={} should in rows={} and ex_col={} should in cols={} " \
                              "with isf_square_dict={} but not, pls check".format(ex_row, isf_square_dict.get("rows"),
                                                                                  ex_col, isf_square_dict.get("cols"),
                                                                                  isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            else:  # output_mode = ""
                return True, isf_square_dict  # bingo(4/4)！

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
    {"data": {"ϵ0": 1, "ϵ1": 1, "ϵ4": -1, "ϵ14": -1,
                  "ϵ5": 1, "ϵ15": -1, "ϵ6": -1, "ϵ16": 1},
     "flags": {"ϵ0": (1, 1), "ϵ1": (1, 1), "ϵ4": (6, 6), "ϵ14": (6, 6),
               "ϵ5": (6, 4), "ϵ15": (3, 1), "ϵ6": (1, 3), "ϵ16": (4, 6)}}
    <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[2, 1, 1]_2.pkl
    {"data": {"ϵ0": 1, "ϵ1": -1, "ϵ4": -1, "ϵ14": 1,
              "ϵ5": -1, "ϵ15": -1, "ϵ6": 1, "ϵ16": 1},
     "flags": {"ϵ0": (1, 4), "ϵ1": (4, 1), "ϵ4": (6, 3), "ϵ14": (3, 6),
               "ϵ5": (6, 1), "ϵ15": (6, 1), "ϵ6": (1, 6), "ϵ16": (1, 6)}}
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


class ΣMDataHelper(object):
    """这里封装供模块内部，σ、μ循环应该取得的数据"""

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
        self.in_matrix_σ_dict = None
        self.in_matrix_μ_dict = None
        self._init_data()

        # self.ϵ_dict_dict = {}
        # self.ϵ_flags_dict = {}
        self.isf_ν_st_list = []
        self.isf_square_dict_dict_by_ν_st = {}
        self.isf_completeness_flag = True  # 默认True，只有需要的时候才会被赋值为False
        self.meta_cgc_ν_τ_list = []

    def __str__(self):
        return "σ={}, μ={} data class".format(self.σ, self.μ)

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
        self.bl_yds_of_σ = self.data_sn_cls.bl_yd_list_dict[tuple(self.σ)]
        self.bl_yds_of_μ = self.data_sn_cls.bl_yd_list_dict[tuple(self.μ)]

        # 3, （tuple(σ'), tuple(μ')）为key， cg_series_st_list为value的dict
        flag, cg_series_st_list_dict = self._load_cg_series_list_dict_by_combination_of_bl_yds()
        if not flag:
            err_msg = "get cg_series_st_list_dict meet error with msg={}".format(cg_series_st_list_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.cg_series_st_list_dict = cg_series_st_list_dict

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

    # def register_ϵ_dict(self, ν, τ, ϵ_dict):
    #     key = tuple(ν) if τ is None else (tuple(ν), τ)
    #     self.ϵ_dict_dict[key] = ϵ_dict
    #
    # def get_ϵ_dict(self, ν, τ):
    #     key = tuple(ν) if τ is None else (tuple(ν), τ)
    #     return self.ϵ_dict_dict.get(key, {})
    #
    # def register_ϵ_flags(self, ν, τ, ϵ_flags):
    #     key = tuple(ν) if τ is None else (tuple(ν), τ)
    #     self.ϵ_flags_dict[key] = ϵ_flags
    #
    # def get_ϵ_flags(self, ν, τ):
    #     key = tuple(ν) if τ is None else (tuple(ν), τ)
    #     return self.ϵ_flags_dict.get(key, {})

    def register_isf_ν_st_list(self, ν_st):
        self.isf_ν_st_list.append(ν_st)

    def get_isf_ν_st_list(self):
        return self.isf_ν_st_list

    def register_isf_square_dict_by_ν_st(self, ν_st, isf_square_dict):
        self.isf_square_dict_dict_by_ν_st[tuple(ν_st)] = isf_square_dict

    def get_isf_square_dict_by_ν_st(self, ν_st):
        return self.isf_square_dict_dict_by_ν_st[tuple(ν_st)]

    def register_meta_isf_completeness_flag_false(self):
        self.isf_completeness_flag = False

    def register_meta_cgc_ν_τ_list(self, ν, τ):
        if (ν, τ) not in self.meta_cgc_ν_τ_list:
            self.meta_cgc_ν_τ_list.append((ν, τ,))

    def get_meta_cgc_ν_τ_list(self):
        return self.meta_cgc_ν_τ_list


class DataHelper(object):
    """这里封装供模块内部，Si循环应该取得的数据"""

    def __init__(self, s_k):
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
        self._init_data()

        # 需手动加载
        self.offset_dict_list = None

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

    def add_more_data(self, yt_st_num_list, yd_st_list):
        # 无法仅通过Sn，还需要传入其他参数才能创建的数据
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


class CalcHelper(object):
    """ISFHelper, CGCHelper公有的函数"""

    def __init__(self):
        self.s_n = None
        self.s_t = None

    def enable_now_s_n(self, s_n):
        if not isinstance(s_n, int) or s_n < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_n={} must be int and >= {}".format(s_n, min(min_s_n_of_isf, min_s_n_of_cgc)))
        self.s_n = s_n
        self.s_t = s_n - 1

    # @staticmethod
    # def _calc_m_with_m_st(yd_st, m_st, bl_yd_list, yt_st_num_dict):
    #     """通过m'计算m，这回，我们可以使用分支律以及yt_num了; 也可以通过m_st=0，计算偏移量"""
    #     # 是core.young_tableaux.py: quickly_calc_young_table_in_decreasing_page_order的特定情况简化
    #     bl_yd_index = bl_yd_list.index(yd_st)  # 必须成立
    #     # 对yd_st的前所有分支的yt_st_num以及本项m_st求和
    #     m = sum([yt_st_num_dict[tuple(bl_yd)] for bl_yd in bl_yd_list[:bl_yd_index]], m_st)
    #     return m


class ISFHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(ISFHelper, self).__init__()

    def calc_meta_isf_dict(self, row_index_tmp_list, ν_st, data_sn, data_st, data_σ_μ):
        """计算元ISF

        对于元σμ组合（无法通过对称性生成的）中的元ISF，使用如下公式计算：
        公式：
        （主公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-192
        （相位公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-195

        正则ISF索引的全部参数为：σ σ' μ μ' ν τ ν' τ'
        表示：|σ σ'> * |μ μ'> 的结果中，|ντ ν'τ'>，的ISF系数平方

        返回flag, meta_isf_square_dict:
        meta_isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                                "cols":[[ν], ([ν], τ), ...],  # 有自由度tuple，无自由度list
                                "isf": isf_square_matrix}  # np.array/sp.Matrix([len(rows), len(cols)])
        """
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
            return True, isf_2_square_dict

        flag, isf_matrix = self._calc_isf_matrix(row_index_tmp_list, ν_st, data_sn, data_σ_μ)
        if not flag:
            err_msg = "get isf_matrix by row_index_tmp_list={}, ν_st={}, data_sn={}, data_σ_μ={} fail " \
                      "with msg={}".format(row_index_tmp_list, ν_st, data_sn, data_σ_μ, isf_matrix)
            logger.error(err_msg)
            return False, err_msg
        try:
            # Return list of triples (eigenval, multiplicity, eigenspace) eigenval升序
            '''
            p.s.
            Matrix([[1, 0, 2], [0, 3, 0], [2, 0, 1]]).eigenvects() return
            [(-1, 1, [Matrix([[-1], [ 0], [ 1]])]),
             (3,  2, [Matrix([[0],  [1],  [0]]),
                      Matrix([[1],  [0],  [1]])])]
            注意这里的结果：单根是正交未归一的；多根不一定正交，也不归一
            '''
            # 到这里isf_matrix是可以通过符号计算，得到无误差的矩阵的。
            # 但后面，求本征值的时候，理论上不一定还能保持符号计算了（一元五次方程无公式解）
            # 不用怕，我们这里还有个终极大招，就是用np数值化计算理论上必然是整数的本征值，再回来使用符号计算
            eigen_tuple_list = isf_matrix.eigenvects()
        except Exception as e:
            logger.error(Exception(e))
            err_msg = "calc eigenvects meet fail with isf_matrix={} s_n={}".format(isf_matrix, self.s_n)
            logger.error(err_msg)
            return False, err_msg

        λ_ν_st = data_st.get_eigenvalue(ν_st)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("@@@@ isf_matrix={} with data_σ_μ={}".format(isf_matrix, data_σ_μ))
            logger.warning("@@@@ row_index_tmp_list={}".format(row_index_tmp_list))
            logger.warning("@@@@ ν_st={}, λ_ν_st={}".format(ν_st, λ_ν_st))
            logger.warning("@@@@ eigen_tuple_list={}".format(eigen_tuple_list))

        # 开始计算ISF
        isf_square_tmp_dict = {"ν_tmp": [],
                               "τ_tmp": [],
                               "isf_tmp": []}
        for e_value, τ_max, e_vectors in eigen_tuple_list[::-1]:  # 按照e_value降序循环
            if τ_max > 1:
                logger.debug("meet multiplicity free cases with e_value={}, τ_max={}, e_vectors={} "
                             "by data_σ_μ={}, ν_st={}".format(e_value, τ_max, e_vectors, data_σ_μ, ν_st))
            # 计算本征值对应的ν
            λ_ν = e_value + λ_ν_st
            ν = self._calc_ν_by_λ_and_bl(λ_ν, data_sn.eigenvalue_list, data_sn.yd_list,
                                         data_sn.bl_yd_list_dict, ν_st)

            # 这里是用[::-1]倒序而不是正序就是为了对齐书中给的多重根自由度选择
            soe_vectors = self._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if τ_max > 1 \
                else [self._calc_orthogonalization_vector(e_vectors[0])]

            flag, τ_tmp_list, isf_phase_vector_list = \
                self._calc_isf_τ_and_phase_vector_list(soe_vectors, ν, ν_st, row_index_tmp_list, τ_max,
                                                       data_sn, data_st, data_σ_μ)
            if not flag:
                err_msg = "_calc_isf_τ_and_phase_vector_list meet error " \
                          "with soe_vectors={}, ν={}, ν_st={}, row_index_tmp_list={}, τ_max={} " \
                          "data_sn={}, data_st={}, data_σ_μ={}, " \
                          "msg={}".format(soe_vectors, ν, ν_st, row_index_tmp_list, τ_max,
                                          data_sn, data_st, data_σ_μ, τ_tmp_list)
                logger.error(err_msg)
                return False, err_msg

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ τ_tmp_list={}".format(τ_tmp_list))
                logger.warning("@@@@ isf_phase_list={}".format(isf_phase_vector_list))

            for τ_tmp in range(1, τ_max + 1):  # 因为按照soe_vectors对应的τ乱序了，所以需要正序
                if τ_max == 1:
                    τ_tmp = None
                τ_tmp_index = τ_tmp_list.index(τ_tmp)
                phase_vector = isf_phase_vector_list[τ_tmp_index]

                # 计算单列的ISF
                isf_square = sp.Matrix([sp.sign(i) * i**2 for i in phase_vector])

                isf_square_tmp_dict["ν_tmp"].append(ν)
                isf_square_tmp_dict["τ_tmp"].append(τ_tmp)  # τ按照代码逻辑，同ν必然升序，且临近
                isf_square_tmp_dict["isf_tmp"].append(isf_square)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("@@@@ isf_square_tmp_dict={}".format(isf_square_tmp_dict))

        # 按照顺序整理矩阵
        row_index_list = [(i[0], i[1]) if i[2] is None else i for i in row_index_tmp_list]
        if len(isf_square_tmp_dict["ν_tmp"]) != len(row_index_list):  # 先检查理论上要求的方阵
            err_msg = "ν_tmp_list={} with len={} and row_index_list={} with len={} must same, pls check".format(
                isf_square_tmp_dict["ν_tmp"], len(isf_square_tmp_dict["ν_tmp"]),
                row_index_list, len(row_index_list))
            logger.error(err_msg)
            return False, err_msg
        meta_isf_square_dict = {"rows": row_index_list,
                           "cols": [],
                           "isf": sp.zeros(len(row_index_list))}
        for ν in data_sn.yd_list:
            while ν in isf_square_tmp_dict["ν_tmp"]:
                # 找数据
                tmp_index = isf_square_tmp_dict["ν_tmp"].index(ν)
                τ = isf_square_tmp_dict["τ_tmp"][tmp_index]
                single_col_index = ν if τ is None else (ν, τ,)
                isf_square = isf_square_tmp_dict["isf_tmp"][tmp_index]
                # 赋值
                # 用上一轮cols的len巧妙确定当前的index
                meta_isf_square_dict["isf"][:, len(meta_isf_square_dict["cols"])] = isf_square
                meta_isf_square_dict["cols"].append(single_col_index)
                # 删tmp数据
                isf_square_tmp_dict["ν_tmp"].pop(tmp_index)
                isf_square_tmp_dict["τ_tmp"].pop(tmp_index)
                isf_square_tmp_dict["isf_tmp"].pop(tmp_index)

        return True, meta_isf_square_dict

    # def complete_same_table_self_symmetry_isf(self, partial_meta_isf_dict, ν_st, data_sn, data_st, data_σ_μ):
    #     """
    #     补完ISF中同表自共轭ev<0的部分（它们会以全列None的形式传进来）
    #     """
    #     # 解析入参
    #     σ_μ_τ_st_list = partial_meta_isf_dict["rows"]
    #     ν_τ_list = partial_meta_isf_dict["cols"]
    #     partial_isf_square = partial_meta_isf_dict["isf"]
    #     # 公共
    #     is_σ_self_conjugate = (data_σ_μ.σ == data_sn.get_tilde(data_σ_μ.σ))
    #     meta_ϵ_st_dict_dict = {}  # 简化版支持用meta_isf_rows的元素拿ϵ_st
    #     for σ_μ_τ_st in σ_μ_τ_st_list:
    #         meta_σ_st, meta_μ_st, τ_st = σ_μ_τ_st if len(σ_μ_τ_st) == 3 else (*σ_μ_τ_st, None)
    #         flag, ϵ_st_dict = load_ϵ(self.s_t, meta_σ_st, meta_μ_st, ν_st, τ_st, is_flag_true_if_not_s_n=False)
    #         if not flag:
    #             err_msg = "load_ϵ fail by self.s_t={}, meta_σ_st={}, meta_μ_st={}, ν_st={}, τ_st={} " \
    #                       "with msg={}".format(self.s_t, meta_σ_st, meta_μ_st, ν_st, τ_st, ϵ_st_dict)
    #             logger.error(err_msg)
    #         meta_ϵ_st_dict_dict[str(σ_μ_τ_st)] = ϵ_st_dict
    #
    #     for col, ν_τ in enumerate(ν_τ_list):  # 按照列循环
    #         meta_ν, τ = ν_τ if isinstance(ν_τ, tuple) else (ν_τ, None)
    #         if data_sn.get_eigenvalue(meta_ν) <= 0:  # 这里巧妙使用ev>0的列做元，去推ev<0的对应列。所以跳过它，执行>0的元
    #             continue
    #         # 用meta列推导对称列
    #         sym_ν = data_sn.get_tilde(meta_ν)
    #         sym_ν_τ = sym_ν if τ is None else (sym_ν, τ)
    #         sym_col = ν_τ_list.index(sym_ν_τ)
    #         meta_ϵ_dict = data_σ_μ.get_ϵ_dict(meta_ν, τ)
    #         for row, σ_μ_τ_st in enumerate(σ_μ_τ_st_list):
    #             meta_σ_st, meta_μ_st, τ_st = σ_μ_τ_st if len(σ_μ_τ_st) == 3 else (*σ_μ_τ_st, None)
    #             if is_σ_self_conjugate:
    #                 sym_σ_st, sym_μ_st = data_st.get_tilde(meta_σ_st), meta_μ_st
    #             else:
    #                 sym_σ_st, sym_μ_st = meta_σ_st, data_st.get_tilde(meta_μ_st)
    #             sym_σ_μ_τ_st = (sym_σ_st, sym_μ_st) if τ_st is None else (sym_σ_st, sym_μ_st, τ_st)
    #             sym_row = σ_μ_τ_st_list.index(sym_σ_μ_τ_st)
    #             # 正式调换None
    #             if partial_isf_square[sym_row, sym_col] is not None:
    #                 err_msg = "sym_row={}, sym_col={} in partial_isf_square={} for same_table should be None".format(
    #                     sym_row, sym_col, partial_isf_square)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #             if is_σ_self_conjugate:
    #                 ϵ = meta_ϵ_dict["ϵ5"] * meta_ϵ_st_dict_dict[str(σ_μ_τ_st)]["ϵ5"]
    #                 butler_σ_σ_st = data_sn.get_butler_phase_factor(data_σ_μ.σ, meta_σ_st)
    #                 butler_ν_ν_st = data_sn.get_butler_phase_factor(meta_ν, ν_st)
    #                 butler_Λ = butler_σ_σ_st * butler_ν_ν_st
    #             else:
    #                 ϵ = meta_ϵ_dict["ϵ6"] * meta_ϵ_st_dict_dict[str(σ_μ_τ_st)]["ϵ6"]
    #                 butler_μ_μ_st = data_sn.get_butler_phase_factor(data_σ_μ.μ, meta_μ_st)
    #                 butler_ν_ν_st = data_sn.get_butler_phase_factor(meta_ν, ν_st)
    #                 butler_Λ = butler_μ_μ_st * butler_ν_ν_st
    #
    #             if _debug_condition(data_σ_μ, ν_st):
    #                 logger.warning("@@@@ ϵ={}, butler_Λ={}, data_σ_μ={}, ν_st={}".format(
    #                     ϵ, butler_Λ, data_σ_μ, ν_st))
    #
    #             partial_isf_square[sym_row, sym_col] = ϵ * butler_Λ * partial_isf_square[row, col]
    #
    #     full_isf_square_dict = {"rows": σ_μ_τ_st_list,
    #                             "cols": ν_τ_list,
    #                             "isf": partial_isf_square}
    #     return True, full_isf_square_dict

    def _calc_isf_τ_and_phase_vector_list(self, soe_vectors, ν, ν_st, row_index_tmp_list, τ_max,
                                          data_sn, data_st, data_σ_μ):
        """ 将经过施密特正交归一化手续后的本征矢量调整为Yamanouchi相位，同时确定非第一分支的τ
        返回结果是按照soe_vectors对应排序的

        对于abs对应不上的，需要调整到与第一分支对齐的soe_vector

        1，将分支律的第一分支首个非零系数调整为正（绝对相位）
        2，非第一分支的，参考第一分支做相对调整（相对相位）
        """
        τ_tmp_list = []  # 这个tmp是相对于完全完成的τ是没有None的，但它int部分的τ就是最终确定的τ
        phase_vector_list = []

        if data_sn.bl_yd_list_dict[tuple(ν)][0] == ν_st:  # ν_st击中ν的第一分支

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ 绝对相位 ν={} 的第一分支是 ν_st={}".format(ν, ν_st))

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
            ν_st_fbl = data_sn.bl_yd_list_dict[tuple(ν)][0]  # ν_st_fbl就是ν_st对应ν的第一分支，它们σ、μ虽然不同，但σ'μ'ν相同

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ 相对相位 ν={} 的第一分支不是 ν_st={} 而是 ν_st_fbl={}".format(ν, ν_st, ν_st_fbl))

            # 三个m的意义分别是：
            # 1 m_ν：               令当前ν_st（非ν的第一分支）的m'取1时，ν的m；
            # 2 m_ν_st_fbl：        令当前ν_st的第一分支ν_st_st的m''取1时，ν_st_fbl（ν的第一分支）的m'
            #                       (ν_st_st是ν_st的第一分支，可以证明它也是ν_st_fbl的分支，但不一定是第一分支)
            # 3,m_ν_by_m_ν_st_fbl： 按照2中方法，确定m_ν_st_fbl后（也意味着确定了ν_st_fbl的杨盘），
            #                       ν对应的m
            m_ν_st = 1
            m_ν = data_sn.quick_calc_m(m_ν_st, ν_st, ν)
            ν_st_st = data_st.bl_yd_list_dict[tuple(ν_st)][0]
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

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ isf_fbl_row_list={}, isf_fbl_col_list={}, isf_fbl_square={} with "
                               "data_σ_μ.σ={}, data_σ_μ.μ={}, ν_st_fbl={}".format(
                    isf_fbl_row_list, isf_fbl_col_list, isf_fbl_square, data_σ_μ.σ, data_σ_μ.μ, ν_st_fbl))

            # 将第一分支(fbl)的τ分配给当前分支，并充分利用结果，把能确定的phase定下来
            '''
            TODO 看用时决定是否优化确定τ的算法
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
                        self._calc_isf_by_known_isf(soe_vector, row_index_tmp_list, ν_st, m_ν_st, m_ν,
                                                    row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl, ν, data_sn, data_σ_μ)
                    if not flag:
                        err_msg = "calc isf_fbl_another meet error " \
                                  "with soe_vector={}, row_index_tmp_list={}, ν_st={}, m_ν_st={}, m_ν={}, " \
                                  "row={}, ν_st_fbl={}, m_ν_st_fbl={}, m_ν_by_m_ν_st_fbl={}, " \
                                  "ν={}, data_sn={}, data_σ_μ={}, " \
                                  "msg={}".format(soe_vector, row_index_tmp_list, ν_st, m_ν_st, m_ν,
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
                            msg = "isf_fbl_another_square_abs={} not in isf_fbl_square_abs_diff_set={}, pls check" \
                                  "with isf_fbl_another={}, row={}, soe_vector={}, data_sn={}, data_st={}, " \
                                  "data_σ_μ={}".format(isf_fbl_another_square_abs, isf_fbl_square_abs_diff_set,
                                                       isf_fbl_another, row, soe_vector, data_sn, data_st, data_σ_μ)
                            logger.debug(msg)
                            logger.debug("will break and use 4-195c to calc")
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

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ τ_candidate={}".format(τ_candidate))
                logger.warning("@@@@ phase_candidate={}".format(phase_candidate))
                logger.warning("@@@@ soe_vectors={}".format(soe_vectors))

            if flag_soe:
                # 到这里，soe_vectors的τ已经完全确定，phase部分确定，下面补充未确定的phase就可以了
                for τ_tmp_set, phase_tmp, soe_vector in zip(τ_candidate, phase_candidate, soe_vectors):
                    # 检查τ
                    if len(τ_tmp_set) != 1:
                        err_msg = "find τ_tmp_set={} not union with data_sn={}, data_st={}, data_σ_μ={}".format(
                            τ_tmp_set, data_sn, data_st, data_σ_μ)
                        logger.error(err_msg)
                        return False, err_msg, None
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

                        if _debug_condition(data_σ_μ, ν_st):
                            logger.warning("@@@@ first_no_0_isf_fbl_row_index={}".format(first_no_0_isf_fbl_row_index))
                            logger.warning("@@@@ first_no_0_isf_fbl_square={}".format(first_no_0_isf_fbl_square))

                        first_no_0_isf_fbl_row = isf_fbl_row_list[first_no_0_isf_fbl_row_index]
                        flag, isf_fbl_another = \
                            self._calc_isf_by_known_isf(soe_vector, row_index_tmp_list, ν_st, m_ν_st, m_ν,
                                                        first_no_0_isf_fbl_row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                        ν, data_sn, data_σ_μ)
                        if not flag:
                            err_msg = "complete isf_fbl_another meet error " \
                                      "with soe_vector={}, row_index_tmp_list={}, ν_st={}, m_ν_st={}, m_ν={}, " \
                                      "first_no_0_isf_fbl_row={}, ν_st_fbl={}, m_ν_st_fbl={}, m_ν_by_m_ν_st_fbl={}, " \
                                      "ν={}, data_sn={}, data_σ_μ={}, " \
                                      "msg={}".format(soe_vector, row_index_tmp_list, ν_st, m_ν_st, m_ν,
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
                                msg = "isf_fbl_another_square_abs={} not eq abs(first_no_0_isf_fbl_square)={}, pls check" \
                                      "with isf_fbl_another={}, first_no_0_isf_fbl_row={}, soe_vector={}, " \
                                      "data_sn={}, data_st={}, data_σ_μ={}".format(
                                    isf_fbl_another_square_abs, first_no_0_isf_fbl_square,
                                    isf_fbl_another, first_no_0_isf_fbl_row, soe_vector, data_sn, data_st, data_σ_μ)
                                logger.debug(msg)
                                logger.debug("will break and use 4-195c to calc")
                                flag_soe = False
                                break
                        # TODO 这里应该是拿待定相位的τ算出来的isf_fbl_another和确定的isf_fbl_another比较啊
                        # phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(soe_vector[first_no_0_isf_fbl_row_index])
                        phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(first_no_0_isf_fbl_square)
                        phase_vector_list.append(phase_tmp_new * soe_vector)

            if not flag_soe:
                τ_tmp_list = []
                phase_vector_list = []

                logger.warning("@@@@ renew isf for ν={}, ν_st={}, data_sn={}, data_st={}, data_σ_μ={}".format(
                    ν, ν_st, data_sn, data_st, data_σ_μ))

                # 使用4-195c重新正算phase_vector，放弃seo_vector
                for τ_tmp in range(1, τ_max + 1):  # 此时，列按照fbl的就是正序了，且它必须有非零τ
                    τ_tmp_list.append(τ_tmp)
                    col = (ν, τ_tmp)
                    isf_fbl_square_by_col = isf_fbl_square[:, isf_fbl_col_list.index(col)]
                    isf_fbl_col_vector = sp.Matrix([sp.sign(i) * sp.sqrt(abs(i)) for i in isf_fbl_square_by_col])
                    phase_vector_tmp = []
                    for row in row_index_tmp_list:
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

                logger.warning("@@@@ renew isf phase_vector_list={}".format(phase_vector_list))
                logger.warning("@@@@ old soe_vectors={}".format(soe_vectors))

                # TODO delete it(因为这部分回推，对计算是多余的，它只是初次写这段代码的检查)
                # 利用new的，再回推fbl，反正以下自洽
                for τ_tmp in range(1, τ_max + 1):
                    col = (ν, τ_tmp)
                    first_isf_fbl_row = isf_fbl_row_list[0]
                    first_isf_fbl_square = isf_fbl_square[:, isf_fbl_col_list.index(col)][0]
                    phase_vector = phase_vector_list[τ_tmp - 1]
                    flag, isf_fbl_another = \
                        self._calc_isf_by_known_isf(phase_vector, row_index_tmp_list, ν_st, m_ν_st, m_ν,
                                                    first_isf_fbl_row, ν_st_fbl, m_ν_st_fbl, m_ν_by_m_ν_st_fbl,
                                                    ν, data_sn, data_σ_μ)
                    if sp.sign(isf_fbl_another) * isf_fbl_another**2 != first_isf_fbl_square:
                        logger.error("$$$$ 回推不成立！")
                        err_msg = "first_isf_fbl_square={}, isf_fbl_another={}".format(
                            first_isf_fbl_square, isf_fbl_another)
                        logger.error(err_msg)
                logger.warning("$$$$ 回推成立！")

        return True, τ_tmp_list, phase_vector_list

    def _calc_isf_by_known_isf(self, k_isf_vector, k_isf_rows, k_ν_st, k_m_ν_st, k_m_ν,
                               u_isf_row, u_ν_st, u_m_ν_st, u_m_ν_by_k_m_ν_st,
                               ν, data_sn, data_σ_μ):
        """
        根据式4-195c，根据已知ISF计算未知

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
        in_matrix_ν_m_m_element = in_matrix_ν[k_m_ν - 1, u_m_ν_by_k_m_ν_st - 1]
        u_isf_element = sum_3_loop / in_matrix_ν_m_m_element
        return True, u_isf_element

    # def _calc_isf_fbl_another_way(self, ν, m_ν, m_ν_by_m_ν_st_fbl, isf_fbl_row, ν_st_fbl,
    #                               row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_σ_μ):
    #     """
    #     isf除了有书中介绍的两种根据St_ISF计算的递推方法外，还可以根据已有Sn的isf，推断一部分isf
    #     见式子 4-195c
    #     这里，是用待定符号的当前isf，计算前面已经确定符号和数值的isf_fbl，其结果最多只可能相差一个正负号
    #     用符号是否相同，就可以判断相对相位是否需要整体改变符号
    #     用数值是否相同，可以判断τ，也可以校验结果的正确性"""
    #     σ_st_fbl, μ_st_fbl = isf_fbl_row[0], isf_fbl_row[1]
    #     τ_st_fbl = isf_fbl_row[2] if len(isf_fbl_row) == 3 else None
    #     flag, cgc_square_st_fbl_dict = load_cgc(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, τ_st_fbl, m_ν_by_m_ν_st_fbl,
    #                                             is_flag_true_if_not_s_n=False)
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("@@@@ cgc_square_st_fbl_dict={}".format(cgc_square_st_fbl_dict))
    #
    #     if not flag:
    #         err_msg = "get cgc_square_dict with self.s_t={}, σ_st_fbl={}, μ_st_fbl={}, " \
    #                   "ν_st_fbl={}, τ_st_fbl={}, cgc_square_st_fbl_dict={} meet error with " \
    #                   "msg={}".format(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, τ_st_fbl,
    #                                   m_ν_by_m_ν_st_fbl, cgc_square_st_fbl_dict)
    #         logger.error(err_msg)
    #         return False, err_msg
    #
    #     sum_3_loop = 0
    #     offset_σ_fbl = data_sn.get_offset(σ_st_fbl, data_σ_μ.σ)
    #     offset_μ_fbl = data_sn.get_offset(μ_st_fbl, data_σ_μ.μ)
    #     for (σ_st, μ_st, τ_st), soe_vector_element in zip(row_index_tmp_list, soe_vector):
    #         flag, cgc_square_st_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st, is_flag_true_if_not_s_n=False)
    #         if not flag:
    #             err_msg = "get cgc_square_st_dict with self.s_t={}, σ_st={}, μ_st={}, ν_st={}, τ_st={}, m_ν_st={} " \
    #                       "meet error with msg={}".format(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st, cgc_square_st_dict)
    #             logger.error(err_msg)
    #             return False, err_msg
    #         sum_3_loop_part = 0
    #         offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
    #         offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
    #         for (m_σ_st_fbl, m_μ_st_fbl), cgc_square_st_fbl_element \
    #                 in cgc_square_st_fbl_dict.items():
    #             m_σ_fbl = offset_σ_fbl + m_σ_st_fbl
    #             m_μ_fbl = offset_μ_fbl + m_μ_st_fbl
    #             for (m_σ_st, m_μ_st), cgc_square_st_element in cgc_square_st_dict.items():
    #                 m_σ = offset_σ + m_σ_st
    #                 m_μ = offset_μ + m_μ_st
    #                 in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_fbl - 1, m_σ - 1]
    #                 in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_fbl - 1, m_μ - 1]
    #                 sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
    #                                    * sp.sign(cgc_square_st_element) * sp.sign(cgc_square_st_fbl_element) \
    #                                    * sp.sqrt(abs(cgc_square_st_element * cgc_square_st_fbl_element))
    #         sum_3_loop += sum_3_loop_part * soe_vector_element
    #
    #     flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in",
    #                                                is_flag_true_if_not_s_n=False)  # (Sn-1, Sn)的对换
    #     if not flag:
    #         err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
    #                   "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
    #         logger.error(err_msg)
    #         return False, err_msg
    #     in_matrix_ν_m_m_element = in_matrix_ν[m_ν - 1, m_ν_by_m_ν_st_fbl - 1]
    #     fbl_isf = sum_3_loop / in_matrix_ν_m_m_element
    #     return True, fbl_isf

    @staticmethod
    def _get_first_no_0_number_from_vector(vector):
        """提取vector中首个误差范围内非0的数字和index"""
        # 现在使用sympy了，不会有误差
        # for number in vector:
        #     if abs(number) < error_value:
        #         continue
        #     return number
        # return None
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
        这里，为了使存在自由度的结果对齐《群表示论的新途径》，使用了一个tricky的技巧，就是取分母之和表示复杂度，返回复杂度最小的一组"""
        soe_1 = sp.GramSchmidt(multi_vectors, True)
        soe_2 = sp.GramSchmidt(multi_vectors[::-1], True)
        soe_1_denominator, soe_2_denominator = 0, 0
        for mat_1, mat_2 in zip(soe_1, soe_2):
            soe_1_denominator += sum(i.as_numer_denom()[1] for i in mat_1)  # 0的'分母'会被取为1，不影响计算啦
            soe_2_denominator += sum(i.as_numer_denom()[1] for i in mat_2)
        rst = soe_1 if soe_1_denominator <= soe_2_denominator else soe_2

        return rst

    @staticmethod
    def _calc_ν_by_λ_and_bl(λ_ν, eigenvalue_list, yd_list, bl_yd_list_dict, ν_st):
        """根据本征值航道，得到ν"""
        count_λ_ν = eigenvalue_list.count(λ_ν)
        if count_λ_ν == 0:
            err_msg = "λ_ν={} should in eigenvalue_list={} but not, pls check".format(λ_ν, eigenvalue_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        elif count_λ_ν == 1:
            return yd_list[eigenvalue_list.index(λ_ν)]

        else:  # 多于1个，需要用"航道"确定
            rst_list = []
            for ev, yd_candidate in zip(eigenvalue_list, yd_list):
                if ev != λ_ν:
                    continue
                bl_yds = bl_yd_list_dict[tuple(yd_candidate)]
                if ν_st in bl_yds:  # 核心就是判断yd_candidate的分支律中是否含有ν_st
                    rst_list.append(yd_candidate)
            if len(rst_list) == 1:
                return rst_list[0]
            else:
                # TODO 这里先报错，暂不确定是否会出现这种情况
                msg = "get rst_list={} len > 1, pls check, " \
                      "with λ_ν={}, eigenvalue_list={}, yd_list={}, bl_yd_list_dict={}, ν_st={}".format(
                    rst_list, λ_ν, eigenvalue_list, yd_list, bl_yd_list_dict, ν_st)
                logger.warning(msg)
                raise Exception(msg)

    @staticmethod
    def calc_row_indexes_tmp(ν_st, data_st, data_σ_μ):
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, τ'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        row_index_tmp_list = []  # [(bl_σ, bl_μ, τ'), ...]
        for bl_σ, bl_μ in product(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ):
            single_ν_cg_series = data_σ_μ.get_cg_series_st(bl_σ, bl_μ, ν_st, data_st.yd_list)
            if single_ν_cg_series == 0:
                continue
            part_rst_list = [(bl_σ, bl_μ, τ_st) if single_ν_cg_series >= 2 else (bl_σ, bl_μ, None)
                             for τ_st in range(1, single_ν_cg_series + 1)]  # [1, 2, ..., single_ν_cg_series]
            row_index_tmp_list += part_rst_list
        return row_index_tmp_list  # 当σ * μ中不包含ν_st时，结果为[]

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

    # def _cgc_st_2_cgc_m_dict(self, cgc_st_tuple, yt_st_num_dict, data_σ_μ):
    #     """取得cg_st_dict，并且，根据分支律，将m_st转化为对应的m
    #     注意：返回的字典，仅仅是把m_st按照分支律准化成了它的m对应，value仍然是st的cgc"""
    #     σ_st, μ_st, ν_st, τ_st, _ = cgc_st_tuple
    #     rst_dict = {}
    #     flag, cgc_st_square_dict = self._load_cgc_with_m1_lru(self.s_t, σ_st, μ_st, ν_st, τ_st)
    #     if not flag:
    #         err_msg = "load_cgc fail by s_t={}, σ_st={}, μ_st={}, ν_st={}, τ_st={} with msg={}".format(
    #             self.s_t, σ_st, μ_st, ν_st, τ_st, cgc_st_square_dict)
    #         logger.error(err_msg)
    #         return False, err_msg
    #     cgc_st_square_n = cgc_st_square_dict.pop("N")
    #     for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
    #         m_σ = self._calc_m_with_m_st(σ_st, m_σ_st, data_σ_μ.bl_yds_of_σ, yt_st_num_dict)
    #         m_μ = self._calc_m_with_m_st(μ_st, m_μ_st, data_σ_μ.bl_yds_of_μ, yt_st_num_dict)
    #         rst_dict[(m_σ, m_μ,)] = cgc_st_square
    #     cgc_st_square_dict["N"] = cgc_st_square_n
    #     rst_dict["N"] = cgc_st_square_n
    #
    #     return True, rst_dict
    
    def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right, data_sn, data_σ_μ):
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

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
                    σ_in_element = data_σ_μ.in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
                    μ_in_element = data_σ_μ.in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1, m_μ_right - 1]
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += in_element_sum * left_tmp * right_tmp

        return True, matrix_element

    # def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right, yt_st_num_dict, data_σ_μ):
    #     """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
    #     matrix_element = 0
    #     # 将St的cgc根据分支律上升到Sn，并整形成py序
    #     # 左
    #     _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_left, yt_st_num_dict, data_σ_μ)
    #     left_n = factor_cgc_left_dict.pop("N")
    #     left_tmp_dict = {}
    #     for left_key, factor_cgc_left in factor_cgc_left_dict.items():
    #         left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
    #         left_tmp_dict[left_key] = left_tmp
    #     factor_cgc_left_dict["N"] = left_n  # 还原字典，避免deepcopy
    #     # 右
    #     if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
    #         right_tmp_dict = left_tmp_dict
    #         left_right_sqrt_n = left_n
    #     else:
    #         _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_right, yt_st_num_dict, data_σ_μ)
    #         right_n = factor_cgc_right_dict.pop("N")
    #         right_tmp_dict = {}
    #         for right_key, factor_cgc_right in factor_cgc_right_dict.items():
    #             right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
    #             right_tmp_dict[right_key] = right_tmp
    #         factor_cgc_right_dict["N"] = right_n
    #         left_right_sqrt_n = sp.sqrt(left_n * right_n)
    #
    #     # 计算matrix_element
    #     for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
    #         for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
    #             in_element_sum = 0
    #             for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
    #                 σ_in_element = data_σ_μ.in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
    #                 μ_in_element = data_σ_μ.in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1, m_μ_right - 1]
    #                 in_element_sum += σ_in_element * μ_in_element
    #             matrix_element += in_element_sum * left_tmp * right_tmp
    #     matrix_element = matrix_element / left_right_sqrt_n
    #
    #     return True, matrix_element

    def _calc_isf_matrix(self, row_index_tmp_list, ν_st, data_sn, data_σ_μ):
        """计算ISF的本征矩阵"""
        # TODO 根据时间反馈决定要不要上多线程/进程
        matrix_div = len(row_index_tmp_list)
        isf_matrix = sp.zeros(matrix_div)

        # 构建主对角元
        for rc, (σ_st, μ_st, τ_st) in enumerate(row_index_tmp_list):
            cgc_st_tuple = (σ_st, μ_st, ν_st, τ_st, 1)
            flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ)
            if not flag:
                err_msg = "get isf_matrix_element by cgc_st_tuple={}, cgc_st_tuple={}, data_sn={}, data_σ_μ={}" \
                          "with msg={}".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                logger.error(err_msg)
                return False, err_msg
            isf_matrix[rc, rc] = isf_matrix_element

        # 构建其他矩阵元素
        if matrix_div >= 2:
            for (row, (σ_st_left, μ_st_left, τ_st_left)), (col, (σ_st_right, μ_st_right, τ_st_right)) \
                    in combinations(enumerate(row_index_tmp_list), 2):
                cgc_st_tuple_left = (σ_st_left, μ_st_left, ν_st, τ_st_left, 1)
                cgc_st_tuple_right = (σ_st_right, μ_st_right, ν_st, τ_st_right, 1)
                flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple_left, cgc_st_tuple_right,
                                                                         data_sn, data_σ_μ)
                if not flag:
                    err_msg = "get isf_matrix_element by cgc_st_tuple_left={}, cgc_st_tuple_right={}, " \
                              "data_sn={}, data_σ_μ={} " \
                              "with msg={}".format(cgc_st_tuple_left, cgc_st_tuple_right, data_sn,
                                                   data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                isf_matrix[row, col] = isf_matrix_element
                isf_matrix[col, row] = isf_matrix_element  # 因为ISF的本征矩阵共轭

        return True, isf_matrix

    # def calc_diff_table_self_symmetry_isf(self, meta_isf_square_dict, meta_ν_st, data_sn, data_st, data_σ_μ):
    #     """
    #     计算异表自共轭ISF
    #     """
    #     if data_σ_μ.σ == data_sn.get_tilde(data_σ_μ.σ):
    #         ϵ_mode = "ϵ5"
    #     else:
    #         ϵ_mode = "ϵ6"
    #     return self.calc_symmetry_isf_include_save(meta_isf_square_dict, meta_ν_st, data_sn, data_st, data_σ_μ,
    #                                                point_ϵ_mode=ϵ_mode)

    def calc_symmetry_isf_include_save(self, meta_isf_square_dict, meta_ν_st, data_sn, data_st, data_σ_μ):
        """
        3，计算σμ外对称ISF
        metaISF循环：(元ISF的标志是ν')
        3.0，数据准备
        3.1，根据ϵ1算¥[μ][σ][ν']¥
        3.2，根据ϵ4算¥[σ~][μ~][ν‘]¥
        3.3，根据ϵ14算¥[μ~][σ~][ν’]¥
        3.4，根据ϵ5算¥[σ~][μ][ν‘~]¥
        3.5，根据ϵ15算¥[μ][σ~][ν’~]¥
        3.6，根据ϵ6算¥[σ][μ~][ν‘~]¥
        3.7，根据ϵ16算¥[μ~][σ][ν’~]¥
        """
        # 公共
        meta_isf_rows = meta_isf_square_dict["rows"]
        meta_isf_cols = meta_isf_square_dict["cols"]
        meta_isf_square_matrix = meta_isf_square_dict["isf"]
        tilde_isf_cols_unsort = [(data_sn.get_tilde(i[0]), i[1]) if isinstance(i, tuple) else data_sn.get_tilde(i)
                                  for i in meta_isf_cols]
        tilde_isf_cols, tilde_meta_2_sym_relationship_cols, _ = self._sort_isf_cols(tilde_isf_cols_unsort, data_sn)
        div = len(meta_isf_rows)  # =len(meta_isf_cols)
        ϵ_st_simple_dict = {}  # 简化版支持用meta_isf_rows的元素拿ϵ_st
        for single_row_st in meta_isf_rows:
            σ_st, μ_st, τ_st = single_row_st if len(single_row_st) == 3 else (*single_row_st, None)
            flag, ϵ_st_dict = load_ϵ(self.s_t, σ_st, μ_st, meta_ν_st, τ_st, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_ϵ fail by self.s_t={}, σ_st={}, μ_st={}, meta_ν_st={}, τ_st={} " \
                          "with msg={}".format(self.s_t, σ_st, μ_st, meta_ν_st, τ_st, ϵ_st_dict)
                logger.error(err_msg)
                return False, err_msg
            ϵ_st_simple_dict[str(single_row_st)] = ϵ_st_dict  # 它们meta_ν_st相同，所以只用σ_st, μ_st, τ_st就能区分
        ϵ_simple_dict = {}  # 简化版支持用meta_isf_cols的元素拿ϵ
        for single_col in meta_isf_cols:
            ν, τ = single_col if isinstance(single_col, tuple) else (single_col, None)
            ϵ_dict = data_σ_μ.get_ϵ_dict(ν, τ)
            ϵ_simple_dict[str(single_col)] = ϵ_dict  # 它们σ, μ相同，所以只用ν, τ就能区分

        # if _debug_condition(data_σ_μ, meta_ν_st):
        #     logger.warning("@@@@ ϵ_st_simple_dict={}, ϵ_simple_dict={}".format(ϵ_st_simple_dict, ϵ_simple_dict))

        # 计算7个对称ISF
        for ϵ_mode in ["ϵ1", "ϵ4", "ϵ14", "ϵ5", "ϵ15", "ϵ6", "ϵ16"]:
            if ϵ_mode == "ϵ1":
                # 3.1，根据ϵ1算¥[μ][σ][ν']¥
                ϵ_tuple = (self.s_n, data_σ_μ.μ, data_σ_μ.σ, meta_ν_st)
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(i[1], i[0]) if len(i) == 2 else (i[1], i[0], i[2]) for i in meta_isf_rows]
            elif ϵ_mode == "ϵ4":
                # 3.2，根据ϵ4算¥[σ~][μ~][ν‘]¥
                ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_sn.get_tilde(data_σ_μ.μ), meta_ν_st)
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(data_st.get_tilde(i[0]), data_st.get_tilde(i[1])) if len(i) == 2 else
                                            (data_st.get_tilde(i[0]), data_st.get_tilde(i[1]), i[2])
                                            for i in meta_isf_rows]
            elif ϵ_mode == "ϵ14":
                # 3.3，根据ϵ14算¥[μ~][σ~][ν’]¥
                ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_sn.get_tilde(data_σ_μ.σ), meta_ν_st)
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(data_st.get_tilde(i[1]), data_st.get_tilde(i[0])) if len(i) == 2 else
                                            (data_st.get_tilde(i[1]), data_st.get_tilde(i[0]), i[2])
                                            for i in meta_isf_rows]
            elif ϵ_mode == "ϵ5":
                # 3.4，根据ϵ5算¥[σ~][μ][ν‘~]¥
                ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_σ_μ.μ, data_st.get_tilde(meta_ν_st))
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(data_st.get_tilde(i[0]), i[1]) if len(i) == 2 else
                                            (data_st.get_tilde(i[0]), i[1], i[2])
                                            for i in meta_isf_rows]
            elif ϵ_mode == "ϵ15":
                # 3.5，根据ϵ15算¥[μ][σ~][ν’~]¥
                ϵ_tuple = (self.s_n, data_σ_μ.μ, data_sn.get_tilde(data_σ_μ.σ), data_st.get_tilde(meta_ν_st))
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(i[1], data_st.get_tilde(i[0])) if len(i) == 2 else
                                            (i[1], data_st.get_tilde(i[0]), i[2])
                                            for i in meta_isf_rows]
            elif ϵ_mode == "ϵ6":
                # 3.6，根据ϵ6算¥[σ][μ~][ν‘~]¥
                ϵ_tuple = (self.s_n, data_σ_μ.σ, data_sn.get_tilde(data_σ_μ.μ), data_st.get_tilde(meta_ν_st))
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(i[0], data_st.get_tilde(i[1])) if len(i) == 2 else
                                            (i[0], data_st.get_tilde(i[1]), i[2])
                                            for i in meta_isf_rows]
            elif ϵ_mode == "ϵ16":
                # 3.7，根据ϵ16算¥[μ~][σ][ν’~]¥
                ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_σ_μ.σ, data_st.get_tilde(meta_ν_st))
                _, is_calc_ed = is_isf_exist(*ϵ_tuple)
                if is_calc_ed is True:
                    continue
                symmetry_isf_rows_unsort = [(data_st.get_tilde(i[1]), i[0]) if len(i) == 2 else
                                            (data_st.get_tilde(i[1]), i[0], i[2])
                                            for i in meta_isf_rows]
            else:
                err_msg = "get unsupported ϵ={}".format(ϵ_mode)
                logger.error(err_msg)
                return False, err_msg

            isf_start_time = time.time()
            symmetry_isf_rows, meta_2_sym_relationship_rows, _ = self._sort_isf_rows(symmetry_isf_rows_unsort,
                                                                                     data_st)
            meta_2_sym_relationship_cols = list(range(div)) if ϵ_mode in ["ϵ1", "ϵ4", "ϵ14"] else \
                tilde_meta_2_sym_relationship_cols

            symmetry_isf_square_matrix = \
                self._calc_symmetry_isf_square_matrix(meta_isf_square_matrix, meta_isf_rows, meta_isf_cols, div,
                                                      ϵ_st_simple_dict, ϵ_simple_dict, ϵ_mode,
                                                      meta_2_sym_relationship_rows, meta_2_sym_relationship_cols,
                                                      meta_ν_st, data_sn, data_σ_μ)

            # symmetry_isf_cols = meta_isf_cols if ϵ_mode in ["ϵ1", "ϵ4", "ϵ14"] else tilde_isf_cols_unsort
            symmetry_isf_cols = meta_isf_cols if ϵ_mode in ["ϵ1", "ϵ4", "ϵ14"] else tilde_isf_cols
            symmetry_isf_square = {"rows": symmetry_isf_rows,
                                   "cols": symmetry_isf_cols,
                                   "isf": symmetry_isf_square_matrix}

            if _debug_condition(data_σ_μ, meta_ν_st):
                logger.warning("@@@@ symmetry_isf_square={} meta_2_sym_relationship_cols={} "
                               "for ϵ_mode={}".format(symmetry_isf_square, meta_2_sym_relationship_cols, ϵ_mode))

            isf_speed_time = int(time.time() - isf_start_time)
            flag, msg = save_isf(*(*ϵ_tuple, symmetry_isf_square, isf_speed_time))
            if not flag:
                err_msg = "save_isf meet error with ϵ_mode={}, s_i={}, σ={}, μ={}, ν_st={}, " \
                          "meta_isf_square_dict={}, " \
                          "msg={}".format(*(ϵ_mode, *ϵ_tuple, meta_isf_square_dict, msg))
                logger.error(err_msg)
                return False, err_msg

        return True, None

    @staticmethod
    def _sort_isf_rows(symmetry_isf_rows_unsort, data_st):
        """将未按照Yamanouchi排序的symmetry_isf_rows_unsort按照Yamanouchi排序
        技巧是利用yd_st_index达到排序yd_st的目的"""
        sort_isf_rows = []
        sym_2_meta_relationship = []  # 形如[1, 3, 2]表示已排序的rows是从meta中哪里来的
        all_yd_st_list = data_st.yd_list
        σ_st_index_list = [all_yd_st_list.index(i[0]) for i in symmetry_isf_rows_unsort]
        σ_st_index_one_list = list(set(σ_st_index_list))
        for σ_st_index_one in sorted(σ_st_index_one_list):  # 排序金牌的σ
            σ_st = all_yd_st_list[σ_st_index_one]
            μ_st_index_list = [all_yd_st_list.index(i[1]) for i in symmetry_isf_rows_unsort if i[0] == σ_st]
            μ_st_index_one_list = list(set(μ_st_index_list))
            for μ_st_index_one in sorted(μ_st_index_one_list):  # 排序银牌的μ
                μ_st = all_yd_st_list[μ_st_index_one]
                # 实有τ的才会是real list，没有的会是[]
                τ_st_list = [i[2] for i in symmetry_isf_rows_unsort if len(i) == 3 and i[0] == σ_st and i[1] == μ_st]
                if τ_st_list:
                    for τ_st in sorted(τ_st_list):  # 排序铜牌的τ
                        single_tuple = (σ_st, μ_st, τ_st)
                        sort_isf_rows.append(single_tuple)
                        sym_2_meta_relationship.append(symmetry_isf_rows_unsort.index(single_tuple))
                else:
                    single_tuple = (σ_st, μ_st)
                    sort_isf_rows.append(single_tuple)
                    sym_2_meta_relationship.append(symmetry_isf_rows_unsort.index(single_tuple))
        # 形如[1, 3, 2]表示meta会到已排序的rows哪里去
        meta_2_sym_relationship = [sym_2_meta_relationship.index(i) for i in range(len(sym_2_meta_relationship))]
        return sort_isf_rows, meta_2_sym_relationship, sym_2_meta_relationship

    @staticmethod
    def _sort_isf_cols(tilde_isf_cols_unsort, data_sn):
        """将未按照Yamanouchi排序的tilde_isf_cols_unsort按照Yamanouchi排序
        技巧是利用yd_index达到排序yd的目的"""
        sort_isf_cols = []
        sym_2_meta_relationship = []  # 形如[1, 3, 2]表示已排序的cols是从meta中哪里来的
        all_yd_list = data_sn.yd_list
        ν_index_list = [all_yd_list.index(i[0]) if isinstance(i, tuple) else all_yd_list.index(i)
                        for i in tilde_isf_cols_unsort]
        ν_index_one_list = list(set(ν_index_list))
        for ν_index_one in sorted(ν_index_one_list):  # 排序金牌的ν
            ν = all_yd_list[ν_index_one]
            # 实有τ的才会是real list，没有的会是[]
            τ_list = [i[1] for i in tilde_isf_cols_unsort if isinstance(i, tuple) and i[0] == ν]
            if τ_list:
                for τ in sorted(τ_list):  # 排序银牌的τ
                    single_tuple = (ν, τ)
                    sort_isf_cols.append(single_tuple)
                    sym_2_meta_relationship.append(tilde_isf_cols_unsort.index(single_tuple))
            else:
                sort_isf_cols.append(ν)
                sym_2_meta_relationship.append(tilde_isf_cols_unsort.index(ν))
        # 形如[1, 3, 2]表示meta会到已排序的cols哪里去
        meta_2_sym_relationship = [sym_2_meta_relationship.index(i) for i in range(len(sym_2_meta_relationship))]
        return sort_isf_cols, meta_2_sym_relationship, sym_2_meta_relationship

    @staticmethod
    def _calc_symmetry_isf_square_matrix(meta_isf_square_matrix, meta_isf_rows, meta_isf_cols, div,
                                         ϵ_st_simple_dict, ϵ_simple_dict, ϵ_key,
                                         meta_2_sym_relationship_rows, meta_2_sym_relationship_cols,
                                         meta_ν_st, data_sn, data_σ_μ):
        symmetry_isf_square_matrix = sp.zeros(div)
        for (row_i, row_st), (col_i, col) in product(enumerate(meta_isf_rows), enumerate(meta_isf_cols)):
            sym_row_i, sym_col_i = meta_2_sym_relationship_rows[row_i], meta_2_sym_relationship_cols[col_i]
            ϵ_st = ϵ_st_simple_dict[str(row_st)][ϵ_key]
            ϵ_sn = ϵ_simple_dict[str(col)][ϵ_key]
            ϵ = ϵ_sn * ϵ_st
            if ϵ_key == "ϵ1":
                # $[μ][σ][ν']$
                butler_Λ = 1
            elif ϵ_key in ["ϵ4", "ϵ14"]:
                # $[σ~][μ~][ν‘]$
                butler_σ_σ_st = data_sn.get_butler_phase_factor(data_σ_μ.σ, row_st[0])
                butler_μ_μ_st = data_sn.get_butler_phase_factor(data_σ_μ.μ, row_st[1])
                butler_Λ = butler_σ_σ_st * butler_μ_μ_st
            elif ϵ_key in ["ϵ5", "ϵ15"]:
                # $[σ~][μ][ν‘~]$
                butler_σ_σ_st = data_sn.get_butler_phase_factor(data_σ_μ.σ, row_st[0])
                butler_ν_ν_st = data_sn.get_butler_phase_factor(col if isinstance(col, list) else col[0], meta_ν_st)
                butler_Λ = butler_σ_σ_st * butler_ν_ν_st
            elif ϵ_key in ["ϵ6", "ϵ16"]:
                # $[σ][μ~][ν‘~]$
                butler_μ_μ_st = data_sn.get_butler_phase_factor(data_σ_μ.μ, row_st[1])
                butler_ν_ν_st = data_sn.get_butler_phase_factor(col if isinstance(col, list) else col[0], meta_ν_st)
                butler_Λ = butler_μ_μ_st * butler_ν_ν_st
            else:
                err_msg = "get unsupported ϵ={}".format(ϵ_key)
                logger.error(err_msg)
                return False, err_msg

            # if _debug_condition(data_σ_μ, meta_ν_st):
            #     logger.warning("@@@@ {}={}, {}_st={}, butler_Λ={}, row_st={}, col={}".format(
            #         ϵ_key, ϵ_sn, ϵ_key, ϵ_st, butler_Λ, row_st, col))

            symmetry_isf_square_matrix[sym_row_i, sym_col_i] = ϵ * butler_Λ * meta_isf_square_matrix[row_i, col_i]

        if _debug_condition(data_σ_μ, meta_ν_st):
            logger.warning("@@@@ symmetry_isf_square_matrix={}, data_σ_μ={}, meta_ν_st={}".format(
                symmetry_isf_square_matrix, data_σ_μ, meta_ν_st))

        return symmetry_isf_square_matrix

    def is_isf_σ_μ_exist(self, σ, μ):
        """只关注存在性，不真正读取数据"""
        if not all(isinstance(yd, list) for yd in [σ, μ]):
            err_msg = "all [σ={}, μ={}] must be list but type [{}, {}]".format(
                σ, μ, type(σ), type(μ))
            logger.error(err_msg)
            return False, err_msg

        flag, file_name = get_isf_file_name(self.s_n, σ, μ, [1])
        if not flag:
            err_msg = "cannot get file_name by s_n={}, σ={}, μ={}, ν_st={} because {}".format(
                self.s_n, σ, μ, [1], file_name)
            logger.error(err_msg)
            return False, err_msg
        folder_name = os.path.dirname(file_name)
        return ISFInfo(self.s_n).folder_exist(folder_name)


class CGCHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(CGCHelper, self).__init__()

    def calc_meta_cgc_include_save(self, meta_isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
        """
        凭借元ISF计算元CGC，并且，计算ϵ
        存储元CGC以及ϵ_dict

        1.2.1，判断，并注册元
        1.2.2，计算元CGC
        """
        # 准备数据
        # 拆解isf_square_dict
        σ_μ_τ_all_st_tuple_list = meta_isf_square_dict["rows"]
        ν_τ_list = meta_isf_square_dict["cols"]
        isf_square_matrix = meta_isf_square_dict["isf"]
        # 得到h_ν_st
        h_ν_st = data_st.get_yt_num(ν_st)

        # 下面的循环，由互不相干的三层构成：行、列、m。更改它们的顺序不影响计算结果。当前使用，先m再列最后行的顺序，比较优
        for m_ν_st in range(1, h_ν_st + 1):
            # 在深层循环前，准备数据，空间换时间
            cgc_st_square_dict_list_by_row = []
            offset_σ_list = []
            offset_μ_list = []
            for σ_μ_τ_all_st in σ_μ_τ_all_st_tuple_list:
                σ_st, μ_st, τ_st = σ_μ_τ_all_st if len(σ_μ_τ_all_st) == 3 else (*σ_μ_τ_all_st, None)
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

            for col, ν_τ in enumerate(ν_τ_list):  # 按照列循环
                # 1.2，根据元ISF计算元CGC【[σ][μ][ν]τ】
                # 拆解列信息
                isf_square_vector = isf_square_matrix.col(col)
                ν, τ = ν_τ if isinstance(ν_τ, tuple) else (ν_τ, None)
                # 每个ν有自己的offset
                offset_ν = data_sn.get_offset(ν_st, ν)
                m_ν = offset_ν + m_ν_st
                # 1.2.1，判断，并注册元信息
                data_σ_μ.register_meta_cgc_ν_τ_list(ν, τ)

                # 1.2.2，已经得到σ,μ,ν,τ,m，开始计算元CGC
                meta_cgc_start_time = time.time()
                meta_cgc_square_dict = {}
                for σ_μ_τ_all_st, isf_square, cgc_st_square_dict, offset_σ, offset_μ in \
                        zip(σ_μ_τ_all_st_tuple_list, isf_square_vector, cgc_st_square_dict_list_by_row,
                            offset_σ_list, offset_μ_list):  # 按照行循环

                    for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():  # 按照CGC_st内容循环
                        cgc_key = (offset_σ + m_σ_st, offset_μ + m_μ_st)

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

                cgc_square_n = sum(abs(i) for i in meta_cgc_square_dict.values())
                if cgc_square_n != 1:  # 归一性检查
                    err_msg = "calc cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, τ={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "meet cgc_square_n={} not eq 1".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν,
                                                                     meta_cgc_square_dict, cgc_square_n)
                    logger.error(err_msg)
                    return False, err_msg

                meta_cgc_speed_time = int(time.time() - meta_cgc_start_time)
                flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, meta_cgc_square_dict,
                                     meta_cgc_speed_time)
                if not flag:
                    err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, τ={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, meta_cgc_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

        return True, None

    # def complete_same_table_self_symmetry_cgc_and_ϵ_include_save(self, complete_isf_square_dict, ν_st,
    #                                                              data_sn, data_st, data_σ_μ):
    #     """
    #     凭借ISF补完元CGC中同表自共轭ev<0的部分，并且，计算那部分的ϵ
    #     存储元CGC以及ϵ_dict
    #     """
    #     # 准备数据
    #     ϵ_speed_time = 0
    #     # 解析入参
    #     ν_τ_list = complete_isf_square_dict["cols"]
    #     # 得到h_ν_st
    #     h_ν_st = data_st.get_yt_num(ν_st)
    #     #
    #     is_σ_self_conjugate = (data_σ_μ.σ == data_sn.get_tilde(data_σ_μ.σ))
    #     is_μ_self_conjugate = (data_σ_μ.μ == data_sn.get_tilde(data_σ_μ.μ))
    #     ϵ_key = "ϵ5" if is_σ_self_conjugate else "ϵ6"
    #
    #     for col, ν_τ in enumerate(ν_τ_list):  # 按照列循环
    #         meta_ν, τ = ν_τ if isinstance(ν_τ, tuple) else (ν_τ, None)
    #         if data_sn.get_eigenvalue(meta_ν) <= 0:  # 这里巧妙使用ev>0的列做元，去推ev<0的对应列。所以跳过它，执行>0的元
    #             continue
    #         sym_ν = data_sn.get_tilde(meta_ν)
    #         h_ν = data_sn.get_yt_num(meta_ν)  # =get_yt_num(sym_ν)
    #         offset_meta_ν = data_sn.get_offset(ν_st, meta_ν)
    #         meta_ϵ_dict = data_σ_μ.get_ϵ_dict(meta_ν, τ)
    #         meta_ϵ = meta_ϵ_dict[ϵ_key]
    #
    #         for m_ν_st in range(1, h_ν_st + 1):
    #             meta_m = offset_meta_ν + m_ν_st
    #             sym_m = h_ν + 1 - meta_m
    #             meta_cgc_tuple = (self.s_n, data_σ_μ.σ, data_σ_μ.μ, meta_ν, τ, meta_m)
    #             sym_cgc_tuple = (self.s_n, data_σ_μ.σ, data_σ_μ.μ, sym_ν, τ, sym_m)
    #             _, is_calc_ed = is_cgc_exist(*sym_cgc_tuple)
    #             if is_calc_ed is True:
    #                 continue
    #
    #             # load元CGC
    #             cgc_start_time = time.time()
    #             flag, meta_cgc_square_dict = load_cgc(*meta_cgc_tuple, is_flag_true_if_not_s_n=False)
    #             if not flag:
    #                 err_msg = "load_cgc fail by meta_cgc_tuple={}, " \
    #                           "with msg={}".format(*meta_cgc_tuple, meta_cgc_square_dict)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #             # 算同表自共轭CGC
    #             if is_σ_self_conjugate:
    #                 sym_cgc_square_dict = \
    #                     self._calc_symmetry_cgc_by_ϵ5(meta_ϵ, meta_cgc_square_dict, meta_ν, meta_m, data_sn, data_σ_μ)
    #             else:
    #                 sym_cgc_square_dict = \
    #                     self._calc_symmetry_cgc_by_ϵ6(meta_ϵ, meta_cgc_square_dict, meta_ν, meta_m, data_sn, data_σ_μ)
    #             # 存
    #             cgc_speed_time = int(time.time() - cgc_start_time)
    #             flag, msg = save_cgc(*(*sym_cgc_tuple, sym_cgc_square_dict, cgc_speed_time))
    #             if not flag:
    #                 err_msg = "save_cgc fail by sym_cgc_tuple={}, sym_cgc_square_dict={} with " \
    #                           "msg={}".format(sym_cgc_tuple, sym_cgc_square_dict, msg)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #             # 计算8个ϵ
    #             if sym_m == h_ν:
    #                 ϵ_start_time_part = time.time()
    #                 sym_ϵ_dict = data_σ_μ.get_ϵ_dict(sym_ν, τ)
    #                 if len(sym_ϵ_dict) == 8:  # ϵ满状态有8个
    #                     continue
    #                 adding_ϵ_dict = self.calc_ϵ_by_m_h(sym_cgc_square_dict, sym_ν, data_sn, data_σ_μ,
    #                                                    is_σ_self_conjugate, is_μ_self_conjugate, False)
    #                 adding_ϵ_dict.update(sym_ϵ_dict)
    #                 data_σ_μ.register_ϵ_dict(sym_ν, τ, adding_ϵ_dict)
    #                 if len(adding_ϵ_dict) == 8:  # ϵ满状态有8个
    #                     flag, msg = save_ϵ(self.s_n, data_σ_μ.σ, data_σ_μ.μ, sym_ν, τ, adding_ϵ_dict)
    #                     if not flag:
    #                         err_msg = "save_ϵ fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
    #                                   "sym_ν={}, τ={}, adding_ϵ_dict={} with " \
    #                                   "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, sym_ν, τ, adding_ϵ_dict, msg)
    #                         logger.error(err_msg)
    #                         return False, err_msg
    #                 ϵ_speed_time_part = int(time.time() - ϵ_start_time_part)
    #                 ϵ_speed_time += ϵ_speed_time_part
    #
    #             if sym_m == 1:
    #                 ϵ_start_time_part = time.time()
    #                 sym_ϵ_dict = data_σ_μ.get_ϵ_dict(sym_ν, τ)
    #                 if len(sym_ϵ_dict) == 8:  # ϵ满状态有8个
    #                     continue
    #                 adding_ϵ_dict = self.calc_ϵ_by_m_1(sym_cgc_square_dict, sym_ν, data_sn, data_σ_μ,
    #                                                    is_σ_self_conjugate, is_μ_self_conjugate, False)
    #                 data_σ_μ.register_ϵ_dict(sym_ν, τ, adding_ϵ_dict)
    #                 if len(adding_ϵ_dict) == 8:  # ϵ满状态有8个
    #                     flag, msg = save_ϵ(self.s_n, data_σ_μ.σ, data_σ_μ.μ, sym_ν, τ, adding_ϵ_dict)
    #                     if not flag:
    #                         err_msg = "save_ϵ fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
    #                                   "sym_ν={}, τ={}, adding_ϵ_dict={} with " \
    #                                   "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, sym_ν, τ, adding_ϵ_dict, msg)
    #                         logger.error(err_msg)
    #                         return False, err_msg
    #                 ϵ_speed_time_part = int(time.time() - ϵ_start_time_part)
    #                 ϵ_speed_time += ϵ_speed_time_part
    #
    #     return True, ϵ_speed_time

    def calc_symmetry_cgc_include_save(self, meta_ν, τ, data_sn, data_σ_μ):
        """
        利用元CGC和ϵ计算全部对称CGC

        metaCGC循环：(元CGC的标志是ν+τ)
        5.0，数据准备
        5.1，根据ϵ1算【[μ][σ][ν]τ】
        5.2，根据ϵ4算【[σ~][μ~][ν]τ】
        5.3，根据ϵ14算【[μ~][σ~][ν]τ】
        5.4，根据ϵ5算【[σ~][μ][ν~]τ】
        5.5，根据ϵ15算【[μ][σ~][ν~]τ】
        5.6，根据ϵ6算【[σ][μ~][ν~]τ】
        5.7，根据ϵ16算【[μ~][σ][ν~]τ】
        """
        h_ν = data_sn.get_yt_num(meta_ν)
        meta_ϵ_dict = data_σ_μ.get_ϵ_dict(meta_ν, τ)
        for m in range(1, h_ν + 1):
            m_tilde = h_ν + 1 - m
            # 5.0，数据准备
            flag, meta_cgc_square_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, meta_ν, τ, m,
                                                  is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, meta_ν={}, τ={}, m={}, " \
                          "with msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, meta_ν, τ, m, meta_cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg

            # 一起算所有CGC太占内存了，所以还是分开算，反正最后的优化方向是这里不算，使用时再现算
            # 5.1，根据ϵ1算【[μ][σ][ν]τ】
            cgc_start_time = time.time()
            cgc1_tuple = (self.s_n, data_σ_μ.μ, data_σ_μ.σ, meta_ν, τ, m)
            _, is_calc_ed = is_cgc_exist(*cgc1_tuple)
            if is_calc_ed is False:  # 能进行到这里的ϵ_x必须非None
                ϵ1 = meta_ϵ_dict["ϵ1"]
                cgc1_square_dict = self.calc_symmetry_cgc_by_ϵ1(ϵ1, meta_cgc_square_dict)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc1_tuple, cgc1_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc1_tuple={}, cgc1_square_dict={} with " \
                              "msg={}".format(cgc1_tuple, cgc1_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.2，根据ϵ4算【[σ~][μ~][ν]τ】
            cgc_start_time = time.time()
            cgc4_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_sn.get_tilde(data_σ_μ.μ), meta_ν, τ, m)
            _, is_calc_ed = is_cgc_exist(*cgc4_tuple)
            if is_calc_ed is False:
                ϵ4 = meta_ϵ_dict["ϵ4"]
                cgc4_square_dict = self._calc_symmetry_cgc_by_ϵ4(ϵ4, meta_cgc_square_dict, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc4_tuple, cgc4_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc4_tuple={}, cgc4_square_dict={} with " \
                              "msg={}".format(cgc4_tuple, cgc4_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.3，根据ϵ14算【[μ~][σ~][ν]τ】
            cgc_start_time = time.time()
            cgc14_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_sn.get_tilde(data_σ_μ.σ), meta_ν, τ, m)
            _, is_calc_ed = is_cgc_exist(*cgc14_tuple)
            if is_calc_ed is False:
                ϵ14 = meta_ϵ_dict["ϵ14"]
                cgc14_square_dict = self._calc_symmetry_cgc_by_ϵ14(ϵ14, meta_cgc_square_dict, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc14_tuple, cgc14_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc14_tuple={}, cgc14_square_dict={} with " \
                              "msg={}".format(cgc14_tuple, cgc14_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.4，根据ϵ5算【[σ~][μ][ν~]τ】
            cgc_start_time = time.time()
            cgc5_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_σ_μ.μ, data_sn.get_tilde(meta_ν), τ, m_tilde)
            _, is_calc_ed = is_cgc_exist(*cgc5_tuple)
            if is_calc_ed is False:
                ϵ5 = meta_ϵ_dict.get("ϵ5")
                cgc5_square_dict = self._calc_symmetry_cgc_by_ϵ5(ϵ5, meta_cgc_square_dict, meta_ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc5_tuple, cgc5_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc5_tuple={}, cgc5_square_dict={} with " \
                              "msg={}".format(cgc5_tuple, cgc5_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.5，根据ϵ15算【[μ][σ~][ν~]τ】
            cgc_start_time = time.time()
            cgc15_tuple = (self.s_n, data_σ_μ.μ, data_sn.get_tilde(data_σ_μ.σ), data_sn.get_tilde(meta_ν), τ, m_tilde)
            _, is_calc_ed = is_cgc_exist(*cgc15_tuple)
            if is_calc_ed is False:
                ϵ15 = meta_ϵ_dict.get("ϵ15")
                cgc15_square_dict = self._calc_symmetry_cgc_by_ϵ15(ϵ15, meta_cgc_square_dict, meta_ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc15_tuple, cgc15_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc15_tuple={}, cgc15_square_dict={} with " \
                              "msg={}".format(cgc15_tuple, cgc15_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.6，根据ϵ6算【[σ][μ~][ν~]τ】
            cgc_start_time = time.time()
            cgc6_tuple = (self.s_n, data_σ_μ.σ, data_sn.get_tilde(data_σ_μ.μ), data_sn.get_tilde(meta_ν), τ, m_tilde)
            _, is_calc_ed = is_cgc_exist(*cgc6_tuple)
            if is_calc_ed is False:
                ϵ6 = meta_ϵ_dict.get("ϵ6")
                cgc6_square_dict = self._calc_symmetry_cgc_by_ϵ6(ϵ6, meta_cgc_square_dict, meta_ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc6_tuple, cgc6_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc6_tuple={}, cgc6_square_dict={} with " \
                              "msg={}".format(cgc6_tuple, cgc6_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 5.7，根据ϵ16算【[μ~][σ][ν~]τ】
            cgc_start_time = time.time()
            cgc16_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_σ_μ.σ, data_sn.get_tilde(meta_ν), τ, m_tilde)
            _, is_calc_ed = is_cgc_exist(*cgc16_tuple)
            if is_calc_ed is False:
                ϵ16 = meta_ϵ_dict.get("ϵ16")
                cgc16_square_dict = self._calc_symmetry_cgc_by_ϵ16(ϵ16, meta_cgc_square_dict, meta_ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*cgc16_tuple, cgc16_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by cgc16_tuple={}, cgc16_square_dict={} with " \
                              "msg={}".format(cgc16_tuple, cgc16_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

        return True, None

    '''
    在以下计算中，需要注意变量的归属问题。
    使用data_σ_μ得到的变量，说明它的计算中，参数是CGC_σμντ的
    而使用ϵ_xxx_tuple的，说明它的计算中，参数是对称目标的
    
    另外，用行列式储存CGC，能得到极大的计算速度和简洁性优化；
    但是，更好的优化是根本不计算它们，只留一个索引，需要调取时，根据索引现场计算，来节省计算和储存资源
    '''
    '''下面的calc_symmetry_cgc_by_ϵ_xxx不怕重复，因为验证过后，它们将不会被计算'''

    @staticmethod
    def calc_symmetry_cgc_by_ϵ1(ϵ1, meta_cgc_square_dict):
        # 根据ϵ1算【[μ][σ][ν]τ m】
        cgc1_square_dict = {(k[1], k[0]): ϵ1 * v for k, v in meta_cgc_square_dict.items()}
        return cgc1_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ4(ϵ4, meta_cgc_square_dict, data_sn, data_σ_μ):
        # 根据ϵ4算【[σ~][μ~][ν]τ m】
        cgc4_square_dict = {}
        h_σ_add_1, h_μ_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(data_σ_μ.μ) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc4_square_dict[(h_σ_add_1 - m_σ, h_μ_add_1 - m_μ)] = ϵ4 * Λ_m_σ * Λ_m_μ * v
        return cgc4_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ14(ϵ14, meta_cgc_square_dict, data_sn, data_σ_μ):
        # 根据ϵ14算【[μ~][σ~][ν]τ m】
        cgc14_square_dict = {}
        h_σ_add_1, h_μ_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(data_σ_μ.μ) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc14_square_dict[(h_μ_add_1 - m_μ, h_σ_add_1 - m_σ)] = ϵ14 * Λ_m_σ * Λ_m_μ * v
        return cgc14_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ5(ϵ5, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ5算【[σ~][μ][ν~]τ m~】
        cgc5_square_dict = {}
        h_σ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            cgc5_square_dict[(h_σ_add_1 - m_σ, m_μ)] = ϵ5 * Λ_m_σ * Λ_m_ν * v
        return cgc5_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ15(ϵ15, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ15算【[μ][σ~][ν~]τ m~】
        cgc15_square_dict = {}
        h_σ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            cgc15_square_dict[(m_μ, h_σ_add_1 - m_σ)] = ϵ15 * Λ_m_σ * Λ_m_ν * v
        return cgc15_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ6(ϵ6, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ6算【[σ][μ~][ν~]τ m~】
        cgc6_square_dict = {}
        h_μ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.μ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc6_square_dict[(m_σ, h_μ_add_1 - m_μ)] = ϵ6 * Λ_m_μ * Λ_m_ν * v
        return cgc6_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ16(ϵ16, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ16算【[μ~][σ][ν~]τ m~】
        cgc16_square_dict = {}
        h_μ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.μ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc16_square_dict[(h_μ_add_1 - m_μ, m_σ)] = ϵ16 * Λ_m_μ * Λ_m_ν * v
        return cgc16_square_dict

    def calc_symmetry_ϵ_by_meta_include_save(self, meta_ν, τ, data_sn, data_σ_μ):
        """
        4，计算σμ外对称ϵ

        从排列组合来看，元ϵ对应的分别是：
        ϵ0    ϵ1    ϵ4    ϵ14   ϵ5    ϵ15   ϵ6    ϵ16
                    _ _   _ _   _       _     _   _
        m1m2  m2m1  m1m2  m2m1  m1m2  m2m1  m1m2  m2m1
        - -   - -                 -   -     -       -
                    Λ1Λ2  Λ2Λ1  Λ1Λν  ΛνΛ1  ΛνΛ2  Λ2Λν
        其中，
        -在下指的是按顺序取小，-在上指的是按顺序取大，Λ对应m的大小和顺序。
        在计算对称的ϵ_dict时，按照公式，实际上拿的也是它们，只不过换了顺序。
        所以，我们只需按照对称性，调整元ϵ_dict中的次序，就可以得到对称ϵ_dict

        4.0，数据准备
        4.1，算【[μ][σ][ν]τ】ϵ1_dict
        4.2，算【[σ~][μ~][ν]τ】ϵ4_dict
        4.3，算【[μ~][σ~][ν]τ】ϵ14_dict
        4.4，算【[σ~][μ][ν~]τ】ϵ5_dict
        4.5，算【[μ][σ~][ν~]τ】ϵ15_dict
        4.6，算【[σ][μ~][ν~]τ】ϵ6_dict
        4.7，算【[μ~][σ][ν~]τ】ϵ16_dict
        """
        # 4.0，数据准备
        h_σ = data_sn.get_yt_num(data_σ_μ.σ)
        h_μ = data_sn.get_yt_num(data_σ_μ.μ)
        Λ_σ_list = data_sn.get_phase_factor_list(data_σ_μ.σ)
        Λ_tilde_σ_list = data_sn.get_phase_factor_list(data_sn.get_tilde(data_σ_μ.σ))
        Λ_μ_list = data_sn.get_phase_factor_list(data_σ_μ.μ)
        Λ_tilde_μ_list = data_sn.get_phase_factor_list(data_sn.get_tilde(data_σ_μ.μ))
        Λ_ν_list = data_sn.get_phase_factor_list(meta_ν)
        Λ_tilde_ν_list = data_sn.get_phase_factor_list(data_sn.get_tilde(meta_ν))

        # 【[σ][μ][ν]τ】
        meta_ϵ_dict = data_σ_μ.get_ϵ_dict(meta_ν, τ)
        meta_ϵ_flags = data_σ_μ.get_ϵ_flags(meta_ν, τ)

        if len(meta_ϵ_dict) != 8 or len(meta_ϵ_flags) != 8:
            err_msg = "meta_ϵ_dict={} and meta_ϵ_flags={} must len8 " \
                      "with meta_ν={}, τ={}".format(meta_ϵ_dict, meta_ϵ_flags, meta_ν, τ)
            logger.error(err_msg)
            return False, err_msg
        # 拿【[σ][μ][ν]τ】验证一下
        sym_ϵ = "ϵ0"
        sym_ϵ_flags = meta_ϵ_flags
        ϵ_map = {"ϵ0": sym_ϵ,
                 "ϵ1": "ϵ1",
                 "ϵ4": "ϵ4",
                 "ϵ14": "ϵ14",
                 "ϵ5": "ϵ5",
                 "ϵ15": "ϵ15",
                 "ϵ6": "ϵ6",
                 "ϵ16": "ϵ16"}
        ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                            Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                            Λ_ν_list, Λ_tilde_ν_list)
        if meta_ϵ_dict != ϵ_dict:
            err_msg = "ϵ0 must get same dict again but not, pls check, " \
                      "with meta_ν={}, τ={}, meta_ϵ_dict={}, ϵ_dict={}".format(meta_ν, τ, meta_ϵ_dict, ϵ_dict)
            logger.error(err_msg)
            return False, err_msg

        # 4.1，算【[μ][σ][ν]τ】ϵ1_dict
        ϵ_tuple = (self.s_n, data_σ_μ.μ, data_σ_μ.σ, meta_ν, τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ1"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ0",
                     "ϵ4": "ϵ14",
                     "ϵ14": "ϵ4",
                     "ϵ5": "ϵ16",
                     "ϵ15": "ϵ6",
                     "ϵ6": "ϵ15",
                     "ϵ16": "ϵ5"}
            sym_ϵ_flags = {k: (meta_ϵ_flags[v][1], meta_ϵ_flags[v][0]) for k, v in ϵ_map.items()}

            if _debug_condition(data_σ_μ, meta_ν):
                logger.warning("@@@@ sym_ϵ={}, ϵ_tuple={}, ϵ_map={}, \n"
                               "meta_ϵ_flags={} \n"
                               "sym_ϵ_flags={}".format(sym_ϵ, ϵ_tuple, ϵ_map, meta_ϵ_flags, sym_ϵ_flags))

            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.2，算【[σ~][μ~][ν]τ】ϵ4_dict
        ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_sn.get_tilde(data_σ_μ.μ), meta_ν, τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ4"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ14",
                     "ϵ4": "ϵ0",
                     "ϵ14": "ϵ1",
                     "ϵ5": "ϵ6",
                     "ϵ15": "ϵ16",
                     "ϵ6": "ϵ5",
                     "ϵ16": "ϵ15"}
            sym_ϵ_flags = {k: (h_σ + 1 - meta_ϵ_flags[v][0], h_μ + 1 - meta_ϵ_flags[v][1]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={} with msg={}".format(ϵ_tuple, ϵ_dict, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.3，算【[μ~][σ~][ν]τ】ϵ14_dict
        ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_sn.get_tilde(data_σ_μ.σ), meta_ν, τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ14"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ4",
                     "ϵ4": "ϵ1",
                     "ϵ14": "ϵ0",
                     "ϵ5": "ϵ15",
                     "ϵ15": "ϵ5",
                     "ϵ6": "ϵ16",
                     "ϵ16": "ϵ6"}
            sym_ϵ_flags = {k: (h_μ + 1 - meta_ϵ_flags[v][1], h_σ + 1 - meta_ϵ_flags[v][0]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.4，算【[σ~][μ][ν~]τ】ϵ5_dict
        ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.σ), data_σ_μ.μ, data_sn.get_tilde(meta_ν), τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ5"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ15",
                     "ϵ4": "ϵ6",
                     "ϵ14": "ϵ16",
                     "ϵ5": "ϵ0",
                     "ϵ15": "ϵ1",
                     "ϵ6": "ϵ4",
                     "ϵ16": "ϵ14"}
            sym_ϵ_flags = {k: (h_σ + 1 - meta_ϵ_flags[v][0], meta_ϵ_flags[v][1]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.5，算【[μ][σ~][ν~]τ】ϵ15_dict
        ϵ_tuple = (self.s_n, data_σ_μ.μ, data_sn.get_tilde(data_σ_μ.σ), data_sn.get_tilde(meta_ν), τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ15"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ5",
                     "ϵ4": "ϵ16",
                     "ϵ14": "ϵ6",
                     "ϵ5": "ϵ1",
                     "ϵ15": "ϵ0",
                     "ϵ6": "ϵ14",
                     "ϵ16": "ϵ4"}
            sym_ϵ_flags = {k: (meta_ϵ_flags[v][1], h_σ + 1 - meta_ϵ_flags[v][0]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.6，算【[σ][μ~][ν~]τ】ϵ6_dict
        ϵ_tuple = (self.s_n, data_σ_μ.σ, data_sn.get_tilde(data_σ_μ.μ), data_sn.get_tilde(meta_ν), τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ6"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ16",
                     "ϵ4": "ϵ5",
                     "ϵ14": "ϵ15",
                     "ϵ5": "ϵ4",
                     "ϵ15": "ϵ14",
                     "ϵ6": "ϵ0",
                     "ϵ16": "ϵ1"}
            sym_ϵ_flags = {k: (meta_ϵ_flags[v][0], h_μ + 1 - meta_ϵ_flags[v][1]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        # 4.7，算【[μ~][σ][ν~]τ】ϵ16_dict
        ϵ_tuple = (self.s_n, data_sn.get_tilde(data_σ_μ.μ), data_σ_μ.σ, data_sn.get_tilde(meta_ν), τ)
        _, is_calc_ed = is_ϵ_exist(*ϵ_tuple)
        if is_calc_ed is False:
            sym_ϵ = "ϵ16"
            ϵ_map = {"ϵ0": sym_ϵ,
                     "ϵ1": "ϵ6",
                     "ϵ4": "ϵ15",
                     "ϵ14": "ϵ5",
                     "ϵ5": "ϵ14",
                     "ϵ15": "ϵ4",
                     "ϵ6": "ϵ1",
                     "ϵ16": "ϵ0"}
            sym_ϵ_flags = {k: (h_μ + 1 - meta_ϵ_flags[v][1], meta_ϵ_flags[v][0]) for k, v in ϵ_map.items()}
            ϵ_dict = self._calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                                                Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list,
                                                Λ_ν_list, Λ_tilde_ν_list)
            flag, msg = save_ϵ(*(*ϵ_tuple, ϵ_dict, sym_ϵ_flags))
            if not flag:
                err_msg = "save_ϵ fail by ϵ_tuple={}, ϵ_dict={}, sym_ϵ_flags={} with msg={}".format(
                    ϵ_tuple, ϵ_dict, sym_ϵ_flags, msg)
                logger.error(err_msg)
                return False, err_msg

        return True, None

    @staticmethod
    def _calc_symmetry_ϵ_dict(sym_ϵ, ϵ_map, meta_ϵ_dict, meta_ϵ_flags, sym_ϵ_flags,
                              Λ_σ_list, Λ_tilde_σ_list, Λ_μ_list, Λ_tilde_μ_list, Λ_ν_list, Λ_tilde_ν_list):
        """
        求ϵ_dict通用的计算部分

        思路：
        ϵ^s{k}  # 求由ϵ{s}造的CGC^s的各种ϵ可以写作ϵ^s{k}
        = sign(CGC^s{k}) * ΛΛ^s{m^s|k}  # 把CGC^s看作元，按照4-122定义给出CGC^s的ϵ0～ϵ16，其中k是对称性种类，也是条件
        = sign(ϵ^0{s} * CGC^0{v} * ΛΛ^0{m^v||s}) * ΛΛ^s{m^s|k}  # CGC^s{k}用CGC^0展开，其中，k^-s=v，s(m^v)=m^k
        = ϵ^0{s} * sign(CGC^0{v}) * ΛΛ^0{m^v||s} * ΛΛ^s{m^s|k}  # 打开sign的括号
        = ϵ^0{s} * (ϵ^0{v} / ΛΛ^0{m^v|v}) * ΛΛ^0{m^v||s} * ΛΛ^s{m^s|k}  # 再按定义用v还原CGC^0{v}
        = ϵ^0{s} * (ϵ^0{v} * ΛΛ^0{m^v|v}) * ΛΛ^0{m^v||s} * ΛΛ^s{m^s|k}  # 对于Λ，*就是/
        = ϵ^0{s} * ϵ^0{v} * ΛΛ^0{m^v|v} * ΛΛ^0{m^v||s} * ΛΛ^s{m^s|k}  # 展开所有的括号
        其中，m的^是序号，｜是定义条件，||是使用条件
        """
        # logger.warning("$$$$ sym_ϵ={}".format(sym_ϵ))
        if sym_ϵ in ["ϵ0", "ϵ6"]:
            Λs_σ = Λ_σ_list
        elif sym_ϵ in ["ϵ1", "ϵ15"]:
            Λs_σ = Λ_μ_list
        elif sym_ϵ in ["ϵ4", "ϵ5"]:
            Λs_σ = Λ_tilde_σ_list
        elif sym_ϵ in ["ϵ14", "ϵ16"]:
            Λs_σ = Λ_tilde_μ_list
        else:
            Λs_σ = None
        if sym_ϵ in ["ϵ0", "ϵ5"]:
            Λs_μ = Λ_μ_list
        elif sym_ϵ in ["ϵ1", "ϵ16"]:
            Λs_μ = Λ_σ_list
        elif sym_ϵ in ["ϵ4", "ϵ6"]:
            Λs_μ = Λ_tilde_μ_list
        elif sym_ϵ in ["ϵ14", "ϵ15"]:
            Λs_μ = Λ_tilde_σ_list
        else:
            Λs_μ = None
        if sym_ϵ in ["ϵ0", "ϵ1", "ϵ4", "ϵ14"]:
            Λs_ν = Λ_ν_list
        elif sym_ϵ in ["ϵ5", "ϵ15", "ϵ6", "ϵ16"]:
            Λs_ν = Λ_tilde_ν_list
        else:
            Λs_ν = None
        # logger.warning("$$$$ Λs_σ={}, Λs_μ={}, Λs_ν={}".format(Λs_σ, Λs_μ, Λs_ν))
        ϵ_dict = {}
        for k, v in ϵ_map.items():
            # ϵ^0{s} 简记为ϵ0s
            ϵ0s = meta_ϵ_dict[sym_ϵ]
            # ϵ^0{v} 简记为ϵ0v
            ϵ0v = meta_ϵ_dict[v]
            # ΛΛ^0{m^v|v} 简记为ΛΛ0_v  # v是0的定义条件，所以用元的Λ；v是0的标号，所以用元的m
            if v in ["ϵ0", "ϵ1"]:
                ΛΛ0_v = 1
            elif v in ["ϵ4", "ϵ14"]:
                ΛΛ0_v = Λ_σ_list[meta_ϵ_flags[v][0] - 1] * Λ_μ_list[meta_ϵ_flags[v][1] - 1]
            elif v in ["ϵ5", "ϵ15"]:
                ΛΛ0_v = Λ_σ_list[meta_ϵ_flags[v][0] - 1] * Λ_ν_list[-1]
            elif v in ["ϵ6", "ϵ16"]:
                ΛΛ0_v = Λ_μ_list[meta_ϵ_flags[v][1] - 1] * Λ_ν_list[-1]
            else:
                ΛΛ0_v = None
            # ΛΛ^0{m^v||s} 简记为ΛΛ0_s  # s是0的使用条件，所以用元的Λ；v是0的标号，所以用元的m
            if sym_ϵ in ["ϵ0", "ϵ1"]:
                ΛΛ0_s = 1
            elif sym_ϵ in ["ϵ4", "ϵ14"]:
                ΛΛ0_s = Λ_σ_list[meta_ϵ_flags[v][0] - 1] * Λ_μ_list[meta_ϵ_flags[v][1] - 1]
            elif sym_ϵ in ["ϵ5", "ϵ15"]:
                ΛΛ0_s = Λ_σ_list[meta_ϵ_flags[v][0] - 1] * Λ_ν_list[-1]
            elif sym_ϵ in ["ϵ6", "ϵ16"]:
                ΛΛ0_s = Λ_μ_list[meta_ϵ_flags[v][1] - 1] * Λ_ν_list[-1]
            else:
                ΛΛ0_s = None
            # ΛΛ^s{m^s|k} 简记为ΛΛs_k  # k是s的定义条件，所以用s的Λ；s是s的标号，所以用s的m，但前面sym_ϵ_flags已经把标号转化为了k
            if k in ["ϵ0", "ϵ1"]:
                ΛΛs_k = 1
            elif k in ["ϵ4", "ϵ14"]:
                ΛΛs_k = Λs_σ[sym_ϵ_flags[k][0] - 1] * Λs_μ[sym_ϵ_flags[k][1] - 1]
            elif k in ["ϵ5", "ϵ15"]:
                ΛΛs_k = Λs_σ[sym_ϵ_flags[k][0] - 1] * Λs_ν[-1]
            elif k in ["ϵ6", "ϵ16"]:
                ΛΛs_k = Λs_μ[sym_ϵ_flags[k][1] - 1] * Λs_ν[-1]
            else:
                ΛΛs_k = None
            ϵ_dict[k] = ϵ0s * ϵ0v * ΛΛ0_v * ΛΛ0_s * ΛΛs_k

        return ϵ_dict

    # def calc_ϵ_and_cgc_dict_by_isf_include_save(self, isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
    #     """
    #     凭借ISF计算CGC，并且，
    #     根据对称性，尽可能多地计算ϵ和对应的CGC；
    #     ϵ_dict和cgc_dict都要储存
    #     """
    #     # 拆解isf_square_dict，准备数据
    #     σ_μ_τ_all_st_tuple_list = isf_square_dict["rows"]
    #     ν_τ_list = isf_square_dict["cols"]
    #     isf_square_matrix = isf_square_dict["isf"]
    #
    #     # TODO 所以一种可能的优化思路是，把m_ν_st循环放在最外面。这样，这个ISF的行和列都可以先一步提取做成字典，空间换时间
    #
    #     for i, ν_τ in enumerate(ν_τ_list):  # 按照列，也就是ν+τ开启循环
    #         isf_square_vector = isf_square_matrix.col(i)
    #         ν, τ = ν_τ if isinstance(ν_τ, tuple) else (ν_τ, None)
    #         h_ν_st = data_st.yt_num_dict[tuple(ν_st)]
    #         offset_of_m_ν = self._calc_m_with_m_st(ν_st, 0, data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)
    #
    #         for m_ν_st in range(1, h_ν_st + 1):
    #             single_cgc_part_start_time = time.time()
    #             m_ν = offset_of_m_ν + m_ν_st
    #
    #             cgc_square_dict = {}
    #             for σ_μ_τ_all_st, isf_square in zip(σ_μ_τ_all_st_tuple_list, isf_square_vector):
    #                 # TODO 这个循环里的load_cgc肯定重复了，它和ν_τ无关；但是，按m全取的话，内存压力太大了
    #                 σ_st, μ_st, τ_st = σ_μ_τ_all_st if len(σ_μ_τ_all_st) == 3 else (*σ_μ_τ_all_st, None)
    #                 flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st,
    #                                                     is_flag_true_if_not_s_n=False)
    #                 if not flag:
    #                     err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, τ_st={}, " \
    #                               "m_ν_st={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, τ_st, m_ν_st,
    #                                                          cgc_st_square_dict)
    #                     logger.error(err_msg)
    #                     return False, err_msg
    #                 offset_of_m_σ = self._calc_m_with_m_st(σ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
    #                                                        data_st.yt_num_dict)
    #                 offset_of_m_μ = self._calc_m_with_m_st(μ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
    #                                                        data_st.yt_num_dict)
    #                 for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
    #                     cgc_key = (offset_of_m_σ + m_σ_st, offset_of_m_μ + m_μ_st)
    #
    #                     if cgc_key not in cgc_square_dict:
    #                         single_cgc_square = isf_square * cgc_st_square
    #                         if single_cgc_square != 0:
    #                             cgc_square_dict[cgc_key] = isf_square * cgc_st_square
    #                     else:
    #                         cgc_new_part = sp.sign(isf_square) * sp.sign(cgc_st_square) \
    #                                        * sp.sqrt(abs(isf_square * cgc_st_square))
    #                         cgc_old_part = sp.sign(cgc_square_dict[cgc_key]) \
    #                                        * sp.sqrt(abs(cgc_square_dict[cgc_key]))
    #                         update_cgc_square = sp.sign(cgc_new_part + cgc_old_part) \
    #                                             * (cgc_new_part + cgc_old_part)**2
    #                         if update_cgc_square != 0:
    #                             cgc_square_dict[cgc_key] = update_cgc_square  # 覆盖
    #                         else:
    #                             cgc_square_dict.pop(cgc_key)
    #
    #             cgc_square_n = sum(abs(i) for i in cgc_square_dict.values())
    #             if cgc_square_n != 1:
    #                 err_msg = "calc cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
    #                           "ν={}, τ={}, m_ν={}, cgc_square_dict={} with " \
    #                           "meet cgc_square_n={} not eq 1".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν,
    #                                                                  cgc_square_dict, cgc_square_n)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #             single_cgc_part_speed_time = int(time.time() - single_cgc_part_start_time)
    #             flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, cgc_square_dict,
    #                                  single_cgc_part_speed_time)
    #             if not flag:
    #                 err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
    #                           "ν={}, τ={}, m_ν={}, cgc_square_dict={} with " \
    #                           "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, m_ν, cgc_square_dict, msg)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #     return True, None


class EHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(EHelper, self).__init__()
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

    @staticmethod
    def _spell_ϵ_key(d3, k4):
        meta = "σμν"
        rst = ""
        for d, k in zip(d3, k4):
            rst += "{}{}".format(meta[d], "~" if k is True else "")
        return rst

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
        ϵ_dict = {}
        ϵ_flags = {}
        # 公共参数
        h_σ = data_sn.get_yt_num(data_σ_μ.σ)
        h_μ = data_sn.get_yt_num(data_σ_μ.μ)
        h_ν = data_sn.get_yt_num(ν)

        # 1，一些可以简化处理的情况可以走快速通道
        if data_σ_μ.σ == data_σ_μ.μ == ν:
            # 这里，除了可以由0，1对称出所有6组σμν外，也不需要遍历所有m的CGC去凑一些dict了
            # 所以，先对m_ν=1、m_ν=h_ν做ϵ，再取对称就可以了
            # 1.1，m_ν=1
            flag, cgc_square_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, 1, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_dict with s_n={}, σ={}, μ={}, ν={}, τ={}, m_ν={} meet error with " \
                          "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, 1, cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_0, m_flags=(None, None, 1))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)

            # 1.2，m_ν=h_ν
            flag, cgc_square_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, h_ν, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_dict with s_n={}, σ={}, μ={}, ν={}, τ={}, h_ν={} meet error with " \
                          "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, τ, h_ν, cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_1, m_flags=(None, None, h_ν))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)

            # 1.3，其他
            for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                # 对0和1采用两种对称（(1,3), (2,3)），实现全部ϵ
                sym_1_d3, sym_1_k4 = (d3[2], d3[1], d3[0]), (k4[2], k4[1], k4[0])
                sym_2_d3, sym_2_k4 = (d3[0], d3[2], d3[1]), (k4[0], k4[2], k4[1])
                sym_1_ϵ_key = self._spell_ϵ_key(sym_1_d3, sym_1_k4)
                sym_2_ϵ_key = self._spell_ϵ_key(sym_2_d3, sym_2_k4)
                ϵ_dict[sym_1_ϵ_key], ϵ_flags[sym_1_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]
                ϵ_dict[sym_2_ϵ_key], ϵ_flags[sym_2_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]

            # 1.4，简化情况可以提前返回
            if len(ϵ_dict) == len(ϵ_flags) == 24:
                return ϵ_dict, ϵ_flags
            else:
                err_msg = "len(ϵ_dict) == len(ϵ_flags) must eq 24 but not " \
                          "with ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags)
                return False, err_msg

        # 2，必须要遍历所有m的CGC的其他情况
        # 按照m_ν循环：
        # m_σ=1、m_μ=1、m_σ=h_σ、m_μ=h_μ的4个字典
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
                                                                   σμν_group=self.σμν_0, m_flags=(None, None, m_ν))
                ϵ_dict.update(adding_ϵ_dict)
                ϵ_flags.update(adding_ϵ_flags)
            if m_ν == h_ν:
                adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(cgc_square_dict, data_σ_μ, ν, data_sn,
                                                                   σμν_group=self.σμν_1, m_flags=(None, None, m_ν))
                ϵ_dict.update(adding_ϵ_dict)
                ϵ_flags.update(adding_ϵ_flags)

        # 2.3，单σ=ν或μ=ν还是有简化算法可以榨取
        if data_σ_μ.σ == ν:
            for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                # 对0和1采用σν对称可以算σμν_4、σμν_5
                if d3 == (0, 1, 2):  # 表示σμν（不看～）
                    sym_d3, sym_k4 = (2, 1, 0), (k4[2], k4[1], k4[0])
                else:  # 这里必然是(1, 0, 2)  # 表示μσν（不看～）
                    sym_d3, sym_k4 = (1, 2, 0), (k4[0], k4[2], k4[1])
                sym_ϵ_key = self._spell_ϵ_key(sym_d3, sym_k4)
                ϵ_dict[sym_ϵ_key], ϵ_flags[sym_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]

        if data_σ_μ.μ == ν:
            for ϵ_key, d3, k4 in chain(self.σμν_0, self.σμν_1):
                # 对0和1采用μν对称可以算σμν_2、σμν_3
                if d3 == (0, 1, 2):  # 表示σμν（不看～）
                    sym_d3, sym_k4 = (0, 2, 1), (k4[0], k4[2], k4[1])
                else:  # 这里必然是(1, 0, 2)  # 表示μσν（不看～）
                    sym_d3, sym_k4 = (2, 0, 1), (k4[2], k4[1], k4[0])
                sym_ϵ_key = self._spell_ϵ_key(sym_d3, sym_k4)
                ϵ_dict[sym_ϵ_key], ϵ_flags[sym_ϵ_key] = ϵ_dict[ϵ_key], ϵ_flags[ϵ_key]

        # 2.4，计算剩余4组ϵ
        if any(i[0] not in ϵ_dict for i in self.σμν_2):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_μ_1_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_2, m_flags=(None, 1, None))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)
        if any(i[0] not in ϵ_dict for i in self.σμν_3):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_μ_h_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_3, m_flags=(None, h_μ, None))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)
        if any(i[0] not in ϵ_dict for i in self.σμν_4):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_σ_1_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_4, m_flags=(1, None, None))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)
        if any(i[0] not in ϵ_dict for i in self.σμν_5):
            adding_ϵ_dict, adding_ϵ_flags = self.calc_adding_ϵ(m_σ_h_dict, data_σ_μ, ν, data_sn,
                                                               σμν_group=self.σμν_5, m_flags=(h_σ, None, None))
            ϵ_dict.update(adding_ϵ_dict)
            ϵ_flags.update(adding_ϵ_flags)

        if len(ϵ_dict) == len(ϵ_flags) == 24:
            return ϵ_dict, ϵ_flags
        else:
            err_msg = "len(ϵ_dict) == len(ϵ_flags) must eq 24 but not " \
                      "with ϵ_dict={}, ϵ_flags={}".format(ϵ_dict, ϵ_flags)
            return False, err_msg

    def calc_adding_ϵ(self, square_dict, data_σ_μ, ν, data_sn, σμν_group=None, m_flags=None):
        """根据传入的6个矩阵，按照定义计算ϵ

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
        Λ_list_list = [data_sn.get_phase_factor_list(yd) for yd in yd_σμν]
        m_set_list = [set(j[i] for j in square_dict.keys()) if m_flags[i] is None else (m_flags[i],) for i in range(3)]

        for ϵ_key, d3, k4 in σμν_group:
            # ϵ_key是字符串形式的键，p.s. μσ~ν~； d3是组合的数字表示，p.s. (1, 0, 2)； k4是组合的共轭bool，p.s.
            if ϵ_key in adding_ϵ_dict:
                continue

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

            # 根据前面计算Λ * Λ
            ΛΛ = 1
            for d, k in zip(d3, k4):
                if k == 1:
                    ΛΛ *= Λ_list_list[d][m_tuple[d] - 1]

            square_dict_key = m_tuple[:-1:] if d3[-1] == default_σμν[-1] else m_tuple  # 真CGC是len2key，拼装的len3key
            adding_ϵ_dict[ϵ_key] = sp.sign(square_dict[square_dict_key]) * ΛΛ
            adding_ϵ_flags[ϵ_key] = m_tuple

            # 一二顺位相等时的简化计算
            if yd_σμν[first_σμν_index] == yd_σμν[second_σμν_index]:
                sym_d3, sym_k4 = (d3[1], d3[0], d3[2]), (k4[1], k4[0], k4[2])
                sym_ϵ_key = self._spell_ϵ_key(sym_d3, sym_k4)
                adding_ϵ_dict[sym_ϵ_key] = adding_ϵ_dict[ϵ_key]
                adding_ϵ_flags[sym_ϵ_key] = adding_ϵ_flags[ϵ_key]

        return adding_ϵ_dict, adding_ϵ_flags

    def calc_symmetry_ϵ_by_meta_include_save(self, meta_ϵ_dict, meta_ϵ_flags, meta_data_σ_μ, meta_ν, τ, meta_data_sn):
        """
        根据元ϵ计算对称的ϵ（坐标系的原点变换）

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
        """
        meta_yd_σμν = (meta_data_σ_μ.σ, meta_data_σ_μ.μ, meta_ν)
        for sym_mode, sym_mode_d3, sym_mode_k4 in \
                chain(self.σμν_0, self.σμν_1, self.σμν_2, self.σμν_3, self.σμν_4, self.σμν_5):
            sym_yd_σμν = (meta_yd_σμν[sym_mode_d3[i]] if sym_mode_k4[i] is False else
                          meta_data_sn.get_tilde(meta_yd_σμν[sym_mode_d3[i]]) for i in range(3))
            sym_ϵ_tuple = (self.s_n, *sym_yd_σμν, τ)
            _, is_calc_ed = is_ϵ_exist(*sym_ϵ_tuple)
            if is_calc_ed is True:
                continue
            sym_ϵ_dict = {}
            sym_ϵ_flags = {}
            '''
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
            '''
        pass

    # @staticmethod
    # def calc_ϵ_by_m_h(cgc_square_dict, ν, data_sn, data_σ_μ):
    #     """根据元CGC的m_h计算ϵ
    #
    #     注意：输入的cgc_square_dict的m应该等于h_ν，它需要在外部确认
    #     注意：书《群表示论的新途径》这段是错的，应该使用英文版《Group Representation Theory for Physicists》中的对应公式
    #     ϵ5：公式4-122b
    #     ϵ6：公式4-122c
    #     ϵ15：自行推导，使ϵ1作用在经过ϵ5对称后的CGC上
    #     ϵ16：自行推导，使ϵ1作用在经过ϵ6对称后的CGC上
    #     """
    #     ϵ_dict = {}
    #     ϵ_flags = {}
    #
    #     # 公共条件
    #     phase_vector_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
    #     phase_vector_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
    #     phase_vector_list_ν = data_sn.get_phase_factor_list(ν)
    #     Λ_h_ν = phase_vector_list_ν[-1]
    #     m_σ_set = set(i[0] for i in cgc_square_dict.keys())
    #     m_μ_set = set(i[1] for i in cgc_square_dict.keys())
    #
    #     # ϵ5条件【[σ~][μ][ν~]τ】：m_ν=h_ν，(m_σ max, m_μ min)，m取（m_σ, m_μ, h_ν）
    #     max_m_σ = max(m_σ_set)
    #     m_μ_second_set = set(i[1] for i in cgc_square_dict.keys() if i[0] == max_m_σ)
    #     min_m_μ = min(m_μ_second_set)
    #     Λ_max_m_σ = phase_vector_list_σ[max_m_σ - 1]  # 物理序要转换成py序
    #     ϵ5_key = (max_m_σ, min_m_μ,)
    #     ϵ5 = sp.sign(cgc_square_dict[ϵ5_key] * Λ_max_m_σ * Λ_h_ν)
    #     ϵ_dict["ϵ5"] = ϵ5
    #     ϵ_flags["ϵ5"] = ϵ5_key
    #
    #     # ϵ15条件【[μ][σ~][ν~]τ】：m_ν=h_ν，(m_μ min, m_σ max)，m取（m_σ, m_μ, h_ν）
    #     min_m_μ = min(m_μ_set)
    #     m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if i[1] == min_m_μ)
    #     max_m_σ = max(m_σ_second_set)
    #     Λ_max_m_σ = phase_vector_list_σ[max_m_σ - 1]  # 物理序要转换成py序
    #     Λ_h_ν = phase_vector_list_ν[-1]
    #     ϵ15_key = (max_m_σ, min_m_μ,)
    #     ϵ15 = sp.sign(cgc_square_dict[ϵ15_key] * Λ_max_m_σ * Λ_h_ν)
    #     ϵ_dict["ϵ15"] = ϵ15
    #     ϵ_flags["ϵ15"] = ϵ15_key
    #
    #     # ϵ6条件【[σ][μ~][ν~]τ】：m_ν=h_ν，(m_σ min, m_μ max)，m取（m_σ, m_μ, h_ν）
    #     min_m_σ = min(m_σ_set)
    #     m_μ_second_set = set(i[1] for i in cgc_square_dict.keys() if i[0] == min_m_σ)
    #     max_m_μ = max(m_μ_second_set)
    #     Λ_max_m_μ = phase_vector_list_μ[max_m_μ - 1]  # 物理序要转换成py序
    #     ϵ6_key = (min_m_σ, max_m_μ,)
    #     ϵ6 = sp.sign(cgc_square_dict[ϵ6_key] * Λ_max_m_μ * Λ_h_ν)
    #     ϵ_dict["ϵ6"] = ϵ6
    #     ϵ_flags["ϵ6"] = ϵ6_key
    #
    #     # ϵ16条件【[μ~][σ][ν~]τ】：m_ν=h_ν，(m_μ max, m_σ min)，m取（m_σ, m_μ, h_ν）
    #     max_m_μ = max(m_μ_set)
    #     m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if i[1] == max_m_μ)
    #     min_m_σ = min(m_σ_second_set)
    #     Λ_max_m_μ = phase_vector_list_μ[max_m_μ - 1]  # 物理序要转换成py序
    #     ϵ16_key = (min_m_σ, max_m_μ,)
    #     ϵ16 = sp.sign(cgc_square_dict[ϵ16_key] * Λ_max_m_μ * Λ_h_ν)
    #     ϵ_dict["ϵ16"] = ϵ16
    #     ϵ_flags["ϵ16"] = ϵ16_key
    #
    #     return ϵ_dict, ϵ_flags

