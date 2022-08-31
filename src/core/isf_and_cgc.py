# -*- coding: utf-8 -*-
"""
this code for creating ISF and CGC
采用对称性方法
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）本段代码主要依据置换群CGC的第二种递推计算方法
# ISF：第四章第19节 公式：4-192（用Sn-1的本征方程计算） 4-189bc（用Sn的CG系数计算） 4-196对称性性质
# CGC：第四章第11、12节
# 整体思路是，先求Sn的ISF，再利用Sn的ISF以及Sn-1的CG，拼出Sn的CG


import copy
import json
import time
# import numpy as np
import sympy as sp
from functools import lru_cache, partial
from itertools import product, combinations, combinations_with_replacement
from conf.cgc_config import default_s_n
from conf.cgc_config import min_s_n_of_isf, min_s_n_of_cgc, min_s_n_of_branching_law, min_s_n_of_ϵ
from core.young_diagrams import load_young_diagrams, calc_young_diagram_dagger
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
    cond_1 = (data_σ_μ.σ == [3, 2] and data_σ_μ.μ == [3, 1, 1] and ν_st == [2, 1, 1])
    cond_2 = (data_σ_μ.σ == [3] and data_σ_μ.μ == [2, 1])
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
    if cond_5_1 or cond_5_2:
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

    # 开启循环计算ISF和CGC(本质上是8个参量σ μ ν β σ’ μ’ ν’ β’的循环，就看怎么写更优了)
    data_si = None
    # data_st = None
    # 按照从小到大的顺序，逐个计算s_i的eigenvalues并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        s_i_isf_speed_time = 0
        logger.debug("Si={}".format(s_i))

        if s_i == 1:
            # 特殊情况，Sn=1时，没有ISF，只有CGC
            σ, μ, ν, β, m = [1], [1], [1], None, 1
            cgc_square_dict = {(1, 1,): 1}
            cgc_speed_time = int(time.time() - s_i_start_time)
            flag, msg = save_cgc(s_i, σ, μ, ν, β, m, cgc_square_dict, cgc_speed_time)
            if not flag:
                err_msg = "save_cgc fail by s_i={}, σ={}, μ={}, ν={}, β={}, m={}, cgc_square_dict={} with " \
                          "msg={}".format(s_i, σ, μ, ν, β, m, cgc_square_dict, msg)
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
        第四版算法如下：
        σμ循环：
        一，第一轮ν'循环（元和ϵ）
        1.1，计算元ISF_σσ'μμ'νν'ββ'
        1.1.1，判断，并注册元CGC（为了2.1.0）
        1.2，根据元ISF计算元CGC【[σ][μ][ν]β】
        1.2.1，判断，并注册元CGC（为了2.1.0）
        1.2.2，计算元CGC
        1.2.3，计算7个ϵ
        注意：循环的是ν'，不是元ISF
        注意：自对称，是特殊情况
        二，第二轮ν'循环（非元，也就是可以被对称的那些）
        2.2，计算对称ISF
        2.1，计算对称CGC
        2.1.0，load元CGC
        2.1.1，根据ϵ1算【[μ][σ][ν]β】
        2.1.2，根据ϵ4算【[σ~][μ~][ν]β】
        2.1.3，根据ϵ14算【[μ~][σ~][ν]β】
        2.1.4，根据ϵ5算【[σ~][μ][ν~]β】
        2.1.5，根据ϵ15算【[μ][σ~][ν~]β】
        2.1.6，根据ϵ6算【[σ][μ~][ν~]β】
        2.1.7，根据ϵ16算【[μ~][σ][ν~]β】
        '''
        for σ, μ in combinations_with_replacement(data_si.yd_list, 2):  # [σ], [μ]双循环，组合不排列
            # σ_μ_start_time = time.time()
            logger.debug("σ={}, μ={}".format(σ, μ))

            # σ, μ循环可以得到的数据
            data_σ_μ = ΣMDataHelper(s_i, σ, μ, data_si, data_st)

            # 第一轮ν'循环，计算元ISF、元CGC和ϵ
            for ν_st in data_st.yd_list:  # ν_st比ν在前

                # 1.1，计算元ISF_σσ'μμ'νν'ββ'
                '''
                所有情况枚举：
                a，σμ=μσ，已经被combinations_with_replacement排除了
                b，is_isf_exist为真，说明本σ、μ组合前，就已经被对称性计算过了（跳过）
                c，(σ=σ~ or μ=μ~) and ν'!=ν'~ and ev(ν') < 0（跳过）
                d，row_index_tmp_list=[]的，表示该组合下不存在ISF（跳过）
                e，(σ=σ~ or μ=μ~) and ν'=ν'~，需检查同表下多β的一致性（正常进行）
                f，其他，在这种顺序下能到这里的，必然是元ISF【注意：元ISF内，也有非元CGC】（正常进行）
                '''
                # b
                _, is_calc_ed = is_isf_exist(s_i, σ, μ, ν_st)
                if is_calc_ed:
                    continue
                # c
                if (data_σ_μ.σ == data_si.get_dagger(data_σ_μ.σ) or data_σ_μ.μ == data_si.get_dagger(data_σ_μ.μ)) \
                        and ν_st != data_st.get_dagger(ν_st) and data_st.get_eigenvalue(ν_st) < 0:
                    continue
                # d
                row_index_tmp_list = isf_func.calc_row_indexes_tmp(ν_st, data_st, data_σ_μ)  # 带着None的rows
                if not row_index_tmp_list:
                    continue
                logger.debug("new meta combination σ={}, μ={}, ν_st={}".format(σ, μ, ν_st))

                single_isf_start_time = time.time()

                # e & f
                data_σ_μ.register_meta_isf_ν_st_list(ν_st)
                '''
                isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                                   "cols":[[ν], ([ν], β), ...],  # 有自由度tuple，无自由度list
                                   "isf": isf_square_matrix}  # sp.Matrix([len(rows), len(cols)])
                '''
                flag, meta_isf_square_dict = isf_func.calc_meta_isf_dict(row_index_tmp_list, ν_st,
                                                                         data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_meta_isf_dict meet error by row_index_tmp_list={}, ν_st={}, data_si={}, " \
                              "data_st={}, data_σ_μ={} with msg={}".format(row_index_tmp_list, ν_st, data_si, data_st,
                                                                           data_σ_μ, meta_isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg

                if _debug_condition(data_σ_μ, ν_st):
                    logger.warning("@@@@ meta_isf_square_dict={} with σ={}, μ={}, ν_st={}".format(
                        meta_isf_square_dict, data_σ_μ.σ, data_σ_μ.μ, ν_st))

                single_isf_speed_time = int(time.time() - single_isf_start_time)
                s_i_isf_speed_time += single_isf_speed_time
                flag, msg = save_isf(s_i, σ, μ, ν_st, meta_isf_square_dict, single_isf_speed_time)
                if not flag:
                    err_msg = "save_isf meet error with s_i={}, σ={}, μ={}, ν_st={}, meta_isf_square_dict={}, " \
                              "msg={}".format(s_i, σ, μ, ν_st, meta_isf_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

                # 1.2，根据元ISF计算元CGC【[σ][μ][ν]β】
                flag, msg = cgc_func.calc_meta_cgc_and_ϵ_include_save(meta_isf_square_dict, ν_st, 
                                                                      data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_meta_cgc_and_ϵ_include_save fail by meta_isf_square_dict={}, ν_st={}, " \
                              "data_si={}, data_st={}, data_σ_μ={} with " \
                              "msg={}".format(meta_isf_square_dict, ν_st, data_si, data_st, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 2.1，计算对称ISF
            meta_isf_ν_st_list = data_σ_μ.get_meta_isf_ν_st_list()
            for meta_ν_st in meta_isf_ν_st_list:
                meta_isf_square_dict = data_σ_μ.get_meta_isf_square_dict_by_ν_st(meta_ν_st)
                flag, msg = isf_func.calc_symmetry_isf_dict_include_save(meta_isf_square_dict, meta_ν_st,
                                                                         data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_symmetry_isf_dict_include_save fail by meta_isf_square_dict={}, meta_ν_st={}, " \
                              "data_sn={}, data_st={}, data_σ_μ={} with " \
                              "msg={}".format(meta_isf_square_dict, meta_ν_st, data_si, data_st, data_σ_μ, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.1，计算对称CGC
            meta_cgc_ν_β_list = data_σ_μ.get_meta_cgc_ν_β_list()
            for (ν, β) in meta_cgc_ν_β_list:
                pass

        # s_i_isf_speed_time  # 它是不算保存，求和出来的
        flag, msg = save_isf_finish_s_n(s_i, s_i_isf_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_isf_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg
        # 它是总时间 - s_i_isf_speed_time 算出来的，也就是除了cgc，还包括了保存ISF的'误差'
        s_i_cgc_speed_time = int(time.time() - s_i_start_time - s_i_isf_speed_time)
        flag, msg = save_cgc_finish_s_n(s_i, s_i_cgc_speed_time, is_check_add_one=True)
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

    isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                       "cols": [[ν], ([ν], β), ...],                   # 有自由度元组，无自由度列表, len(rows)=len(cols)
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
        # 输出单行，额外输入行指标([σ'], [μ'], β') or ([σ'], [μ'])
        if not isinstance(ex_params, tuple):
            err_msg = "ex_params={} with type={} must be tuple".format(ex_params, type(ex_params))
            logger.error(err_msg)
            return False, err_msg
    elif output_mode == "single_col":
        # 输出单列，额外输入列指标[ν] or ([ν], β)
        if not isinstance(ex_params, (list, tuple)):
            err_msg = "ex_params={} with type={} must be tuple or list".format(ex_params, type(ex_params))
            logger.error(err_msg)
            return False, err_msg
    elif output_mode == "single_isf":
        # 输出单独isf，额外输入行列指标([σ'], [μ'], β') or ([σ'], [μ']) / [ν] or ([ν], β)
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
                # 输出单行，额外输入行指标([σ'], [μ'], β') or ([σ'], [μ'])
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
                # 输出单列，额外输入列指标[ν] or ([ν], β)
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


def save_cgc(s_n: int, σ: list, μ: list, ν: list, β: (int, None), m: int, cgc_square_data, speed_time: int=0):
    """
    这个db用来存CGC

    <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_β_m.pkl  # 无多重性，则为[ν]_m.pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_β_m,  # 无多重性，则为[ν]_m
    "data": cgc_square_dict,
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    β对应[ν]的多重性;
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
    if β is not None and not isinstance(β, int):
        err_msg = "β={} with type={} must be None or int".format(β, type(β))
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
    _, file_name = get_cgc_file_name(s_n, σ, μ, ν, β, m)

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


def load_cgc(s_n: int, σ: list, μ: list, ν: list, β: (int, None), m: int, is_flag_true_if_not_s_n=True):
    """
    取得s_n下指定[σ][μ][ν] β m 的cg系数字典
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
    if β is not None and (not isinstance(β, int) or β <= 0):
        err_msg = "type(β={})={} must be None or real int".format(β, type(β))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(m, int) or m <= 0:
        err_msg = "type(m={})={} must be real int".format(m, type(m))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, β, m)
    if not flag:
        err_msg = "cannot get file_name by s_n={} σ={}, μ={}, ν={}, β={}, m={}, because {}".format(
            s_n, σ, μ, ν, β, m, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = CGCInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query cgc with s_n={} σ={}, μ={}, ν={}, β={}, m={}, because {}".format(
            s_n, σ, μ, ν, β, m, data)
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


def is_cgc_exist(s_n: int, σ: list, μ: list, ν: list, β: (int, None), m: int):
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
    if β is not None and (not isinstance(β, int) or β <= 0):
        err_msg = "type(β={})={} must be None or real int".format(β, type(β))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(m, int) or m <= 0:
        err_msg = "type(m={})={} must be real int".format(m, type(m))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, β, m)
    # flag, file_name = get_cgc_file_name(s_n, σ, μ, ν, β, m, is_full_path=True)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, σ={}, μ={}, ν={}, β={}, m={} because {}".format(
            s_n, σ, μ, ν, β, m, file_name)
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


def save_ϵ(s_n: int, σ: list, μ: list, ν: list, β: (int, None), ϵ_dict: dict):
    """
    这个db用来存ϵ

    <CG>/ϵ_info/Sn/[σ]_[μ]/[ν]_β.pkl  # 无多重性，则为[ν].pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_β,  # 无多重性，则为[ν]
    "data": ϵ_dict,
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    β对应[ν]的多重性;

    例如，
    <CG>/ϵ_info/S5/[3, 2]_[3, 1, 1]/[3, 1, 1]_2.pkl
    {"ϵ1": int, "ϵ4": int, "ϵ14": int, "ϵ5": int, "ϵ15": int, "ϵ6": int, "ϵ16": int}
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
    if β is not None and not isinstance(β, int):
        err_msg = "β={} with type={} must be None or int".format(β, type(β))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(ϵ_dict, dict):
        err_msg = "ϵ_dict={} with type={} must be dict".format(ϵ_dict, type(ϵ_dict))
        logger.error(err_msg)
        return False, err_msg

    db_info = EInfo(s_n)
    _, file_name = get_ϵ_file_name(s_n, σ, μ, ν, β)

    table = {"file_name": file_name,
             "data": ϵ_dict,
             "flags": {}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


# def update_ϵ(s_n: int, σ: list, μ: list, ν: list, β: (int, None), ϵ_dict: dict):
#     if not isinstance(s_n, int) or s_n <= 0:
#         err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
#         logger.error(err_msg)
#         return False, err_msg
#     if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
#         err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
#             σ, μ, ν, type(σ), type(μ), type(ν))
#         logger.error(err_msg)
#         return False, err_msg
#     if β is not None and not isinstance(β, int):
#         err_msg = "β={} with type={} must be None or int".format(β, type(β))
#         logger.error(err_msg)
#         return False, err_msg
#     if not isinstance(ϵ_dict, dict):
#         err_msg = "ϵ_dict={} with type={} must be dict".format(ϵ_dict, type(ϵ_dict))
#         logger.error(err_msg)
#         return False, err_msg
#
#     db_info = EInfo(s_n)
#     _, file_name = get_ϵ_file_name(s_n, σ, μ, ν, β)
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


def load_ϵ(s_n: int, σ: list, μ: list, ν: list, β: (int, None), is_flag_true_if_not_s_n=True):
    """
    取得s_n下指定[σ][μ][ν] β 的ϵ字典
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
    if β is not None and (not isinstance(β, int) or β <= 0):
        err_msg = "type(β={})={} must be None or real int".format(β, type(β))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_ϵ_file_name(s_n, σ, μ, ν, β)
    if not flag:
        err_msg = "cannot get file_name by s_n={} σ={}, μ={}, ν={}, β={}, because {}".format(
            s_n, σ, μ, ν, β, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = EInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query ϵ with s_n={} σ={}, μ={}, ν={}, β={}, because {}".format(
            s_n, σ, μ, ν, β, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        ϵ_dict = data.get("data")
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

        self.ϵ_dict_dict = {}
        self.meta_isf_ν_st_list = []
        self.meta_isf_square_dict_dict_by_ν_st = {}
        self.meta_cgc_ν_β_list = []

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

    def register_ϵ_dict(self, ν, β, ϵ_dict):
        key = tuple(ν) if β is None else (tuple(ν), β)
        self.ϵ_dict_dict[key] = ϵ_dict

    def get_ϵ_dict(self, ν, β):
        key = tuple(ν) if β is None else (tuple(ν), β)
        return self.ϵ_dict_dict.get(key, {})

    def register_meta_isf_ν_st_list(self, ν_st):
        self.meta_isf_ν_st_list.append(ν_st)

    def get_meta_isf_ν_st_list(self):
        return self.meta_isf_ν_st_list

    def register_meta_isf_square_dict_by_ν_st(self, ν_st, meta_isf_square_dict):
        self.meta_isf_square_dict_dict_by_ν_st[tuple(ν_st)] = meta_isf_square_dict

    def get_meta_isf_square_dict_by_ν_st(self, ν_st):
        return self.meta_isf_square_dict_dict_by_ν_st[tuple(ν_st)]

    def register_meta_cgc_ν_β_list(self, ν, β):
        if (ν, β) not in self.meta_cgc_ν_β_list:
            self.meta_cgc_ν_β_list.append((ν, β,))

    def get_meta_cgc_ν_β_list(self):
        return self.meta_cgc_ν_β_list


class DataHelper(object):
    """这里封装供模块内部，Si循环应该取得的数据"""

    def __init__(self, s_k):
        if not isinstance(s_k, int) or s_k < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_k={} must be int and >= 1".format(s_k))
        self.s_n = s_k  # 既可以实例化Sn，也可以实例化St

        self.yd_list = None
        self.bl_yd_list_dict = None
        self.yt_num_list = None
        self.eigenvalue_list = None
        self.yd_dagger_list = None
        self.phase_factor_list_list = None
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
        for yd in yd_list:
            flag, bl_dict = load_branching_law(self.s_n, yd, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get bl_dict meet error with s_n={}, yd={}, msg={}".format(self.s_n, yd, bl_dict)
                logger.error(err_msg)
                raise Exception(err_msg)
            bl_yd_list_dict[tuple(yd)] = bl_dict.get("before_YD")
        self.bl_yd_list_dict = bl_yd_list_dict

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
        yd_dagger_list = []
        for yd in yd_list:
            flag, yd_dagger = calc_young_diagram_dagger(yd, is_check_yd=False)
            if not flag:
                err_msg = "calc yd_dagger meet error with s_n={}, yd={}, msg={}".format(self.s_n, yd, yd_dagger)
                logger.error(err_msg)
                raise Exception(err_msg)
            yd_dagger_list.append(yd_dagger)
        self.yd_dagger_list = yd_dagger_list

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

    def get_dagger(self, yd):
        return self.yd_dagger_list[self.yd_list.index(yd)]

    def get_phase_factor_list(self, yd):
        return self.phase_factor_list_list[self.yd_list.index(yd)]
    # TODO 查一遍还在用index调取的方法，统一改为函数

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
        """计算ISF

        公式：
        （主公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-192
        （相位公式）《群表示论的新途径》陈金全（上海科学技术出版社1984） 4-195

        正则ISF索引的全部参数为：σ σ' μ μ' ν β ν' β'
        表示：|σ σ'> * |μ μ'> 的结果中，|νβ ν'β'>，的ISF系数平方

        返回flag, isf_square_dict:
        isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                           "cols":[[ν], ([ν], β), ...],  # 有自由度tuple，无自由度list
                           "isf": isf_square_matrix}  # np.array/sp.Matrix([len(rows), len(cols)])
        """
        # TODO σ=σ～ and ν=ν～ and β！=None 时的自检查
        # 只象征性检查Sn、St
        if not all(self.s_n == s_n for s_n in [data_sn.s_n, data_σ_μ.s_n]) \
                or not all(self.s_t == s_t for s_t in [data_st.s_n, data_σ_μ.s_t]):
            err_msg = "input wrong data with s_n={} should eq all(data_sn.s_n={}, data_σ_μ.s_n={}) " \
                      "and s_t={} should eq all(data_st.s_n={}, data_σ_μ.s_t={}) but not".format(
                self.s_n, data_sn.s_n, data_σ_μ.s_n, self.s_t, data_st.s_n, data_σ_μ.s_t)
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
            注意，这里的结果，但根是正交未归一的，多根不一定正交，也不归一
            '''
            # TODO 到这里isf_matrix是可以通过符号计算，得到无误差的矩阵的。
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
            logger.warning("@@@@ isf_matrix={} \n with σ={} μ={}".format(isf_matrix, data_σ_μ.σ, data_σ_μ.μ))
            logger.warning("@@@@ ν_st={}, λ_ν_st={}".format(ν_st, λ_ν_st))
            logger.warning("@@@@ eigen_tuple_list={}".format(eigen_tuple_list))

        # 开始计算ISF
        isf_square_tmp_dict = {"ν_tmp": [],
                               "β_tmp": [],
                               "isf_tmp": []}
        for e_value, β_max, e_vectors in eigen_tuple_list[::-1]:  # 按照e_value降序循环
            # 计算本征值对应的ν
            λ_ν = e_value + λ_ν_st
            ν = self._calc_ν_by_λ_and_bl(λ_ν, data_sn.eigenvalue_list, data_sn.yd_list,
                                         data_sn.bl_yd_list_dict, ν_st)

            # 这里是用[::-1]倒序而不是正序就是为了对齐书中给的多重根自由度选择
            soe_vectors = self._calc_schmidt_orthogonalization_tricky(e_vectors)[::-1] if β_max > 1 \
                else [self._calc_orthogonalization_vector(e_vectors[0])]

            flag, β_tmp_list, isf_phase_list = \
                self._calc_isf_β_and_phase_list(soe_vectors, ν, ν_st, row_index_tmp_list, β_max,
                                                data_sn, data_st, data_σ_μ)
            if not flag:
                err_msg = "calc _calc_isf_β_and_phase_list meet error " \
                          "with soe_vectors={}, ν={}, ν_st={}, row_index_tmp_list={}, β_max={} " \
                          "data_sn={}, data_st={}, data_σ_μ={}, " \
                          "msg={}".format(soe_vectors, ν, ν_st, row_index_tmp_list, β_max, data_sn, data_st, data_σ_μ,
                                          β_tmp_list)
                logger.error(err_msg)
                return False, err_msg

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ β_tmp_list={}".format(β_tmp_list))
                logger.warning("@@@@ isf_phase_list={}".format(isf_phase_list))

            for β_tmp in range(1, β_max + 1):  # 因为按照soe_vectors对应的β乱序了，所以需要正序
                if β_max == 1:
                    β_tmp = None
                β_tmp_index = β_tmp_list.index(β_tmp)
                isf_phase = isf_phase_list[β_tmp_index]
                soe_vector = soe_vectors[β_tmp_index]

                phase_vector = soe_vector * isf_phase

                # 计算单列的ISF
                isf_square = sp.Matrix([sp.sign(i) * i**2 for i in phase_vector])

                isf_square_tmp_dict["ν_tmp"].append(ν)
                isf_square_tmp_dict["β_tmp"].append(β_tmp)  # β按照代码逻辑，同ν必然升序，且临近
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
        isf_square_dict = {"rows": row_index_list,
                           "cols": [],
                           "isf": sp.zeros(len(row_index_list))}
        for ν in data_sn.yd_list:
            while ν in isf_square_tmp_dict["ν_tmp"]:
                # 找数据
                tmp_index = isf_square_tmp_dict["ν_tmp"].index(ν)
                β = isf_square_tmp_dict["β_tmp"][tmp_index]
                single_col_index = ν if β is None else (ν, β,)
                isf_square = isf_square_tmp_dict["isf_tmp"][tmp_index]
                # 赋值
                isf_square_dict["isf"][:, len(isf_square_dict["cols"])] = isf_square  # 用上一轮cols的len巧妙确定当前的index
                isf_square_dict["cols"].append(single_col_index)
                # 删tmp数据
                isf_square_tmp_dict["ν_tmp"].pop(tmp_index)
                isf_square_tmp_dict["β_tmp"].pop(tmp_index)
                isf_square_tmp_dict["isf_tmp"].pop(tmp_index)

        return True, isf_square_dict

    def _calc_isf_β_and_phase_list(self, soe_vectors, ν, ν_st, row_index_tmp_list, β_max,
                                   data_sn, data_st, data_σ_μ):
        """ 将经过施密特正交归一化手续后的本征矢量调整为Yamanouchi相位，同时确定非第一分支的β
        返回结果是按照soe_vectors对应排序的

        1，将分支律的第一分支首个非零系数调整为正（绝对相位）
        2，非第一分支的，参考第一分支做相对调整（相对相位）
        """
        β_tmp_list = []  # 这个tmp是相对于完全完成的β是没有None的，但它int部分的β就是最终确定的β
        phase_list = []

        if data_sn.bl_yd_list_dict[tuple(ν)][0] == ν_st:  # ν_st击中ν的第一分支

            if _debug_condition(data_σ_μ, ν_st):
                logger.warning("@@@@ 绝对相位 ν={} 的第一分支是 ν_st={}".format(ν, ν_st))

            # 当ν_st为ν第一分支的时候，令首个非零系数符号为+，称为绝对相位
            for β_tmp, soe_vector in zip(range(1, β_max + 1), soe_vectors):
                β_tmp = None if β_max == 1 else β_tmp
                _, no_0_element = self._get_first_no_0_number_from_vector(soe_vector)
                if no_0_element is None:
                    err_msg = "find schmidt_orthogonalization_vector={} all approximately equal to 0, pls check".format(
                        soe_vector)
                    logger.error(err_msg)
                    return False, err_msg, None
                β_tmp_list.append(β_tmp)
                absolute_phase = sp.sign(no_0_element)
                phase_list.append(absolute_phase)

            return True, β_tmp_list, phase_list

        else:  # 未击中情况，则相对相位需要参考它的分支律第一分支（ν_st_fbl）的绝对相位
            # 被参考的第一分支被称为fbl: first branching law
            '''
            这里的计算思路是，根据式4-195c，既可以用第一分支，计算其他分支的isf；也可以用其他分支，反向计算第一分支
            所以，如果仅仅差一个相位的isf已经得到，就可以用它反推第一分支，相同则不改相位；不同则*-1
            另外，非第一分支的β也可以用此方法确定（也是反推，直到确定属于哪个β）
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

            # 将第一分支(fbl)的β分配给当前分支，并充分利用结果，把能确定的phase定下来
            '''
            TODO 看用时决定是否优化确定β的算法
            为了区分β，我们采取一个最笨的办法，就是：
            1，维护一个flag_uncertain_β_num，用来判断是否完成目标
            2，维护β_by_soe_vectors，用来表示当前每个vector的β范围（循环结束后应该是len1的）
            3，维护phase_by_soe_vectors，用来充分利用结果，把已经计算出的phase写进去，后面就不用再算了
            for 循环所有soe_vectors的行（fbl作为答案，行号和soe_vectors的是一致的）
                1，判断flag_uncertain_β_num
                2，把当前未确定β的列 对应的isf_fbl_square提取出来
                3，计算abs集合，判断是否需要跳过该结果
                
                for 循环soe_vectors的所有列（未完全确定β的列，它们有的有可能已经被缩小范围了）
                    0，跳过已经确定β的列
                    1，用整体符号待定的isf计算对应行的isf_fbl
                    2，缩小范围
                    3，判断是否全部β都已经被区分
                        3-1，区分，退出
                        3-2，未区分，开始下一轮循环
            对未确定phase的列，补充phase
            '''
            β_fbl_list = [None] if β_max == 1 else list(range(1, β_max + 1))  # load的isf必然按照β升序排列
            # 初始所有β_fbl_list都是候选，随着逐渐缩小范围直至全部确定  # 它与soe_vectors顺序是对应的
            β_candidate = [set(β_fbl_list) for i in range(β_max)]  # 它用列表解析就是深拷贝，用*就是浅拷贝
            phase_candidate = [None] * β_max  # 确定一个填一个，允许循环结束后有未确定的
            fbl_cols_useful = [ν] if β_max == 1 else list((ν, i) for i in β_fbl_list)  # 指的是fbl中能用于判断β和phase的那部分
            fbl_cols_useful_index = [isf_fbl_col_list.index(col) for col in fbl_cols_useful]
            β_make_sure_set = {None} if β_max == 1 else set()  # 也作为统一旗标，避免指标太多，逻辑混乱

            for row_index, row in enumerate(isf_fbl_row_list):  # 既是fbl的，也是soe_vectors的
                if len(β_make_sure_set) == β_max:
                    break
                isf_fbl_square_by_row = isf_fbl_square[row_index, :]
                isf_fbl_square_useful_list = [isf_fbl_square_by_row[i] for i in fbl_cols_useful_index]
                isf_fbl_square_abs_β_set_dict = {}  # abs(isf_fbl_square): set(β_possible)
                for single_fbl_β, single_isf_fbl_square in zip(β_fbl_list, isf_fbl_square_useful_list):
                    if single_fbl_β in β_make_sure_set:
                        continue
                    isf_fbl_square_abs = abs(single_isf_fbl_square)
                    if isf_fbl_square_abs in isf_fbl_square_abs_β_set_dict:
                        isf_fbl_square_abs_β_set_dict[isf_fbl_square_abs] &= {single_fbl_β}
                    else:
                        isf_fbl_square_abs_β_set_dict[isf_fbl_square_abs] = {single_fbl_β}
                isf_fbl_square_abs_diff_set = set(isf_fbl_square_abs_β_set_dict.keys())
                if len(isf_fbl_square_abs_diff_set) == 1:
                    # 绝对值全相等，说明无区分β的能力，跳过
                    # TODO 如果全部跳过，说明β任意，程序暂时未考虑出现此情况
                    # TODO 还有一种未考虑情况如：a1=[1/2, -1/2] a2=[-1/2, 1/2], a3=[1/3, -2/3], a4=[2/3, -1/3]...
                    continue

                for single_β_candidate_set, (col_index, soe_vector) in zip(β_candidate, enumerate(soe_vectors)):
                    if len(single_β_candidate_set) == 1:  # 跳过确定好的β
                        continue
                    # 用当前待定β的列soe_vector，去计算已知的isf_fbl_square
                    flag, isf_fbl_another = self._calc_isf_fbl_another_way(
                        ν, m_ν, m_ν_by_m_ν_st_fbl, row, ν_st_fbl,
                        row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_σ_μ)
                    if not flag:
                        err_msg = "calc _calc_isf_fbl_another_way meet error with ν={}, m_ν={}, " \
                                  "m_ν_by_m_ν_st_fbl={}, row={}, ν_st_fbl={}, row_index_tmp_list={}, " \
                                  "soe_vector={}, ν_st={}, m_ν_st={}, data_sn={}, data_σ_μ={}, " \
                                  "msg={}".format(ν, m_ν, m_ν_by_m_ν_st_fbl, row, ν_st_fbl, row_index_tmp_list,
                                                  soe_vector, ν_st, m_ν_st, data_sn, data_σ_μ, isf_fbl_another)
                        logger.error(err_msg)
                        return False, err_msg, None
                    isf_fbl_another_square_abs = isf_fbl_another ** 2
                    # abs检查，要求必须和某一个已知结果相等
                    # TODO 暂时关掉检查1，统计一些次数
                    if isf_fbl_another_square_abs not in isf_fbl_square_abs_diff_set:
                        err_msg = "$$$$ isf_fbl_another_square_abs={} not in isf_fbl_square_abs_diff_set={}, pls check" \
                                  "with isf_fbl_another={}, row={}, soe_vector={}, data_sn={}, data_st={}, " \
                                  "data_σ_μ={}".format(isf_fbl_another_square_abs, isf_fbl_square_abs_diff_set,
                                                       isf_fbl_another, row, soe_vector, data_sn, data_st, data_σ_μ)
                        logger.warning(err_msg)
                        continue
                        # return False, err_msg, None
                    # 缩小β范围
                    isf_fbl_square_abs_β_set = isf_fbl_square_abs_β_set_dict[isf_fbl_another_square_abs]
                    single_β_candidate_set &= isf_fbl_square_abs_β_set  # 用旧范围交集新范围
                    if len(single_β_candidate_set) == 1:
                        # 完全确定了一个β，flag_uncertain_β_num和isf_fbl_square_abs_β_set_dict都要随之改变
                        isf_fbl_square_abs_β_set_dict[isf_fbl_another_square_abs] -= single_β_candidate_set
                        β_make_sure_set |= single_β_candidate_set
                    # 顺带确定phase（只有完全确定了β的才能定phase）
                    if isf_fbl_another_square_abs != 0 and len(single_β_candidate_set) == 1:
                        # phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(soe_vector[row_index])
                        fbl_col = ν if β_max == 1 else (ν, list(single_β_candidate_set)[0])
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
                logger.warning("@@@@ β_candidate={}".format(β_candidate))
                logger.warning("@@@@ phase_candidate={}".format(phase_candidate))
                logger.warning("@@@@ soe_vectors={}".format(soe_vectors))

            # 到这里，soe_vectors的β已经完全确定，phase部分确定，下面补充未确定的phase就可以了
            for β_tmp_set, phase_tmp, soe_vector in zip(β_candidate, phase_candidate, soe_vectors):
                # 检查β
                if len(β_tmp_set) != 1:
                    err_msg = "find β_tmp_set={} not union with data_sn={}, data_st={}, data_σ_μ={}".format(
                        β_tmp_set, data_sn, data_st, data_σ_μ)
                    logger.error(err_msg)
                    return False, err_msg, None
                β_tmp = list(β_tmp_set)[0]
                β_tmp_list.append(β_tmp)
                # 补漏phase
                if phase_tmp is not None:
                    phase_list.append(phase_tmp)
                else:
                    fbl_col = ν if β_max == 1 else (ν, β_tmp)
                    fbl_col_index = isf_fbl_col_list.index(fbl_col)
                    isf_fbl_square_same_β = isf_fbl_square[:, fbl_col_index]
                    first_no_0_isf_fbl_row_index, first_no_0_isf_fbl_square = \
                        self._get_first_no_0_number_from_vector(isf_fbl_square_same_β)

                    if _debug_condition(data_σ_μ, ν_st):
                        logger.warning("@@@@ first_no_0_isf_fbl_row_index={}".format(first_no_0_isf_fbl_row_index))
                        logger.warning("@@@@ first_no_0_isf_fbl_square={}".format(first_no_0_isf_fbl_square))

                    first_no_0_isf_fbl_row = isf_fbl_row_list[first_no_0_isf_fbl_row_index]
                    flag, isf_fbl_another = self._calc_isf_fbl_another_way(
                        ν, m_ν, m_ν_by_m_ν_st_fbl, first_no_0_isf_fbl_row, ν_st_fbl,
                        row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_σ_μ)
                    if not flag:
                        err_msg = "calc _calc_isf_fbl_another_way meet error with ν={}, m_ν={}, " \
                                  "m_ν_by_m_ν_st_fbl={}, first_no_0_isf_fbl_row={}, ν_st_fbl={}, " \
                                  "row_index_tmp_list={}, soe_vector={}, ν_st={}, m_ν_st={}, " \
                                  "data_sn={}, data_σ_μ={}, " \
                                  "msg={}".format(ν, m_ν, m_ν_by_m_ν_st_fbl, first_no_0_isf_fbl_row,
                                                  ν_st_fbl, row_index_tmp_list, soe_vector, ν_st, m_ν_st,
                                                  data_sn, data_σ_μ, isf_fbl_another)
                        logger.error(err_msg)
                        return False, err_msg, None
                    isf_fbl_another_square_abs = isf_fbl_another ** 2
                    # abs检查
                    # TODO 暂时关掉检查2，统计一些次数
                    if isf_fbl_another_square_abs != abs(first_no_0_isf_fbl_square):
                        err_msg = "$$$$ isf_fbl_another_square_abs={} not eq abs(first_no_0_isf_fbl_square)={}, pls check" \
                                  "with isf_fbl_another={}, first_no_0_isf_fbl_row={}, soe_vector={}, " \
                                  "data_sn={}, data_st={}, data_σ_μ={}".format(
                            isf_fbl_another_square_abs, first_no_0_isf_fbl_square,
                            isf_fbl_another, first_no_0_isf_fbl_row, soe_vector, data_sn, data_st, data_σ_μ)
                        logger.warning(err_msg)
                    # if isf_fbl_another_square_abs != abs(first_no_0_isf_fbl_square):
                    #     err_msg = "isf_fbl_another_square_abs={} not eq abs(first_no_0_isf_fbl_square)={}, pls check" \
                    #               "with isf_fbl_another={}, first_no_0_isf_fbl_row={}, soe_vector={}, " \
                    #               "data_sn={}, data_st={}, data_σ_μ={}".format(
                    #         isf_fbl_another_square_abs, first_no_0_isf_fbl_square,
                    #         isf_fbl_another, first_no_0_isf_fbl_row, soe_vector, data_sn, data_st, data_σ_μ)
                    #     logger.error(err_msg)
                    #     return False, err_msg, None
                    # TODO 这里应该是拿待定相位的β算出来的isf_fbl_another和确定的isf_fbl_another比较啊
                    # phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(soe_vector[first_no_0_isf_fbl_row_index])
                    phase_tmp_new = sp.sign(isf_fbl_another) * sp.sign(first_no_0_isf_fbl_square)
                    phase_list.append(phase_tmp_new)

            return True, β_tmp_list, phase_list

    def _calc_isf_fbl_another_way(self, ν, m_ν, m_ν_by_m_ν_st_fbl, isf_fbl_row, ν_st_fbl,
                                  row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_σ_μ):
        """
        isf除了有书中介绍的两种根据St_ISF计算的递推方法外，还可以根据已有Sn的isf，推断一部分isf
        见式子 4-195c
        这里，是用待定符号的当前isf，计算前面已经确定符号和数值的isf_fbl，其结果最多只可能相差一个正负号
        用符号是否相同，就可以判断相对相位是否需要整体改变符号
        用数值是否相同，可以判断β，也可以校验结果的正确性"""
        σ_st_fbl, μ_st_fbl = isf_fbl_row[0], isf_fbl_row[1]
        β_st_fbl = isf_fbl_row[2] if len(isf_fbl_row) == 3 else None
        flag, cgc_square_st_fbl_dict = load_cgc(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, β_st_fbl, m_ν_by_m_ν_st_fbl,
                                                is_flag_true_if_not_s_n=False)

        if _debug_condition(data_σ_μ, ν_st):
            logger.warning("@@@@ cgc_square_st_fbl_dict={}".format(cgc_square_st_fbl_dict))

        if not flag:
            err_msg = "get cgc_square_dict with self.s_t={}, σ_st_fbl={}, μ_st_fbl={}, " \
                      "ν_st_fbl={}, β_st_fbl={}, cgc_square_st_fbl_dict={} meet error with " \
                      "msg={}".format(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, β_st_fbl,
                                      m_ν_by_m_ν_st_fbl, cgc_square_st_fbl_dict)
            logger.error(err_msg)
            return False, err_msg

        sum_3_loop = 0
        offset_σ_fbl = data_sn.get_offset(σ_st_fbl, data_σ_μ.σ)
        offset_μ_fbl = data_sn.get_offset(μ_st_fbl, data_σ_μ.μ)
        for (σ_st, μ_st, β_st), soe_vector_element in zip(row_index_tmp_list, soe_vector):
            flag, cgc_square_st_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_st_dict with self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, m_ν_st={} " \
                          "meet error with msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, cgc_square_st_dict)
                logger.error(err_msg)
                return False, err_msg
            sum_3_loop_part = 0
            offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
            offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
            for (m_σ_st_fbl, m_μ_st_fbl), cgc_square_st_fbl_element \
                    in cgc_square_st_fbl_dict.items():
                m_σ_fbl = offset_σ_fbl + m_σ_st_fbl
                m_μ_fbl = offset_μ_fbl + m_μ_st_fbl
                for (m_σ_st, m_μ_st), cgc_square_st_element in cgc_square_st_dict.items():
                    m_σ = offset_σ + m_σ_st
                    m_μ = offset_μ + m_μ_st
                    in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_fbl - 1, m_σ - 1]
                    in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_fbl - 1, m_μ - 1]
                    sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
                                       * sp.sign(cgc_square_st_element) * sp.sign(cgc_square_st_fbl_element) \
                                       * sp.sqrt(abs(cgc_square_st_element * cgc_square_st_fbl_element))
            sum_3_loop += sum_3_loop_part * soe_vector_element

        flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in",
                                                   is_flag_true_if_not_s_n=False)  # (Sn-1, Sn)的对换
        if not flag:
            err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
                      "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
            logger.error(err_msg)
            return False, err_msg
        in_matrix_ν_m_m_element = in_matrix_ν[m_ν - 1, m_ν_by_m_ν_st_fbl - 1]
        fbl_isf = sum_3_loop / in_matrix_ν_m_m_element
        return True, fbl_isf

    # def _calc_isf_fbl_another_way(self, ν, m_ν, m_ν_by_m_ν_st_fbl, isf_fbl_row, ν_st_fbl,
    #                               row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ):
    #     """
    #     isf除了有书中介绍的两种根据St_ISF计算的递推方法外，还可以根据已有Sn的isf，推断一部分isf
    #     见式子 4-195c
    #     这里，是用待定符号的当前isf，计算前面已经确定符号和数值的isf_fbl，其结果最多只可能相差一个正负号
    #     用符号是否相同，就可以判断相对相位是否需要整体改变符号
    #     用数值是否相同，可以判断β，也可以校验结果的正确性"""
    #     σ_st_fbl, μ_st_fbl = isf_fbl_row[0], isf_fbl_row[1]
    #     β_st_fbl = isf_fbl_row[2] if len(isf_fbl_row) == 3 else None
    #     flag, cgc_square_st_fbl_dict = load_cgc(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, β_st_fbl, m_ν_by_m_ν_st_fbl,
    #                                             is_flag_true_if_not_s_n=False)
    #
    #     if _debug_condition(data_σ_μ, ν_st):
    #         logger.warning("@@@@ cgc_square_st_fbl_dict={}".format(cgc_square_st_fbl_dict))
    #
    #     if not flag:
    #         err_msg = "get cgc_square_dict with self.s_t={}, σ_st_fbl={}, μ_st_fbl={}, " \
    #                   "ν_st_fbl={}, β_st_fbl={}, cgc_square_st_fbl_dict={} meet error with " \
    #                   "msg={}".format(self.s_t, σ_st_fbl, μ_st_fbl, ν_st_fbl, β_st_fbl,
    #                                   m_ν_by_m_ν_st_fbl, cgc_square_st_fbl_dict)
    #         logger.error(err_msg)
    #         return False, err_msg
    #     cgc_square_fbl_n = cgc_square_st_fbl_dict.pop("N")
    #
    #     sum_3_loop = 0
    #     for (σ_st, μ_st, β_st), soe_vector_element in zip(row_index_tmp_list, soe_vector):
    #         flag, cgc_square_st_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, is_flag_true_if_not_s_n=False)
    #         if not flag:
    #             err_msg = "get cgc_square_st_dict with self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, m_ν_st={} " \
    #                       "meet error with msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, cgc_square_st_dict)
    #             logger.error(err_msg)
    #             return False, err_msg
    #         cgc_square_dict_n = cgc_square_st_dict.pop("N")
    #         sum_3_loop_part = 0
    #         for (m_σ_st_fbl, m_μ_st_fbl), cgc_square_st_fbl_element \
    #                 in cgc_square_st_fbl_dict.items():
    #             m_σ_fbl = self._calc_m_with_m_st(σ_st_fbl, m_σ_st_fbl, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
    #                                              data_st.yt_num_dict)
    #             m_μ_fbl = self._calc_m_with_m_st(μ_st_fbl, m_μ_st_fbl, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
    #                                              data_st.yt_num_dict)
    #             for (m_σ_st, m_μ_st), cgc_square_st_element in cgc_square_st_dict.items():
    #                 m_σ = self._calc_m_with_m_st(σ_st, m_σ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
    #                                              data_st.yt_num_dict)
    #                 m_μ = self._calc_m_with_m_st(μ_st, m_μ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
    #                                              data_st.yt_num_dict)
    #                 in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_fbl - 1, m_σ - 1]
    #                 in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_fbl - 1, m_μ - 1]
    #                 sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
    #                                    * sp.sign(cgc_square_st_element) * sp.sign(cgc_square_st_fbl_element) \
    #                                    * sp.sqrt(abs(cgc_square_st_element * cgc_square_st_fbl_element
    #                                                  / (cgc_square_dict_n * cgc_square_fbl_n)))
    #         sum_3_loop += sum_3_loop_part * soe_vector_element
    #
    #     flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in")  # (Sn-1, Sn)的对换
    #     if not flag:
    #         err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
    #                   "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
    #         logger.error(err_msg)
    #         return False, err_msg
    #     in_matrix_ν_m_m_element = in_matrix_ν[m_ν - 1, m_ν_by_m_ν_st_fbl - 1]
    #     fbl_isf = sum_3_loop / in_matrix_ν_m_m_element
    #     return True, fbl_isf

    # def _calc_fbl_isf_tmp(self):
    #     """TODO 看看有关fblerence这一段能不能独立出来成为一个函数"""
    #     pass

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
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, β'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        row_index_tmp_list = []  # [(bl_σ, bl_μ, β'), ...]
        for bl_σ, bl_μ in product(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ):
            single_ν_cg_series = data_σ_μ.get_cg_series_st(bl_σ, bl_μ, ν_st, data_st.yd_list)
            if single_ν_cg_series == 0:
                continue
            part_rst_list = [(bl_σ, bl_μ, β_st) if single_ν_cg_series >= 2 else (bl_σ, bl_μ, None)
                             for β_st in range(1, single_ν_cg_series + 1)]  # [1, 2, ..., single_ν_cg_series]
            row_index_tmp_list += part_rst_list
        return row_index_tmp_list  # 当σ * μ中不包含ν_st时，结果为[]

    @lru_cache()
    def _load_cgc_with_m1_by_input_json(self, s_n, j_σ, j_μ, j_ν, β):
        """lru版本的load_cgc，通过将入参更改为json.dumps实现序列化/json.loads还原"""
        σ, μ, ν = map(json.loads, [j_σ, j_μ, j_ν])
        return load_cgc_with_m1(s_n, σ, μ, ν, β)

    def _load_cgc_with_m1_lru(self, s_n, σ, μ, ν, β):
        j_σ, j_μ, j_ν = map(json.dumps, [σ, μ, ν])
        return self._load_cgc_with_m1_by_input_json(s_n, j_σ, j_μ, j_ν, β)

    def _cgc_st_2_cgc_m_dict(self, cgc_st_tuple, data_sn, data_σ_μ):
        """取得cg_st_dict，并且，根据分支律，将m_st转化为对应的m
        注意：返回的字典，仅仅是把m_st按照分支律准化成了它的m对应，value仍然是st的cgc"""
        σ_st, μ_st, ν_st, β_st, _ = cgc_st_tuple
        rst_dict = {}
        flag, cgc_st_square_dict = self._load_cgc_with_m1_lru(self.s_t, σ_st, μ_st, ν_st, β_st)
        if not flag:
            err_msg = "load_cgc fail by s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={} with msg={}".format(
                self.s_t, σ_st, μ_st, ν_st, β_st, cgc_st_square_dict)
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
    #     σ_st, μ_st, ν_st, β_st, _ = cgc_st_tuple
    #     rst_dict = {}
    #     flag, cgc_st_square_dict = self._load_cgc_with_m1_lru(self.s_t, σ_st, μ_st, ν_st, β_st)
    #     if not flag:
    #         err_msg = "load_cgc fail by s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={} with msg={}".format(
    #             self.s_t, σ_st, μ_st, ν_st, β_st, cgc_st_square_dict)
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
        for rc, (σ_st, μ_st, β_st) in enumerate(row_index_tmp_list):
            cgc_st_tuple = (σ_st, μ_st, ν_st, β_st, 1)
            flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ)
            if not flag:
                err_msg = "get isf_matrix_element by cgc_st_tuple={}, cgc_st_tuple={}, data_sn={}, data_σ_μ={}" \
                          "with msg={}".format(cgc_st_tuple, cgc_st_tuple, data_sn, data_σ_μ, isf_matrix_element)
                logger.error(err_msg)
                return False, err_msg
            isf_matrix[rc, rc] = isf_matrix_element

        # 构建其他矩阵元素
        if matrix_div >= 2:
            for (row, (σ_st_left, μ_st_left, β_st_left)), (col, (σ_st_right, μ_st_right, β_st_right)) \
                    in combinations(enumerate(row_index_tmp_list), 2):
                cgc_st_tuple_left = (σ_st_left, μ_st_left, ν_st, β_st_left, 1)
                cgc_st_tuple_right = (σ_st_right, μ_st_right, ν_st, β_st_right, 1)
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

    def calc_symmetry_isf_dict_include_save(self, meta_isf_square_dict, meta_ν_st, data_sn, data_st, data_σ_μ):
        """
        利用元ISF和ϵ计算全部对称ISF

        2.1.1，根据ϵ1算\\[μ][σ][ν']\\
        2.1.2，根据ϵ4算\\[σ~][μ~][ν‘]\\
        2.1.3，根据ϵ14算\\[μ~][σ~][ν’]\\
        2.1.4，根据ϵ5算\\[σ~][μ][ν‘~]\\
        2.1.5，根据ϵ15算\\[μ][σ~][ν’~]\\
        2.1.6，根据ϵ6算\\[σ][μ~][ν‘~]\\
        2.1.7，根据ϵ16算\\[μ~][σ][ν’~]\\
        """
        # 公共
        meta_isf_rows = meta_isf_square_dict["rows"]
        meta_isf_cols = meta_isf_square_dict["cols"]
        meta_isf_square_matrix = meta_isf_square_dict["isf"]
        dagger_isf_cols_unsort = [(data_st.get_dagger(i[0]), i[1]) if isinstance(i, tuple) else data_sn.get_dagger(i)
                                  for i in meta_isf_cols]
        dagger_isf_cols, dagger_meta_2_sym_relationship_cols, _ = self._sort_isf_cols(dagger_isf_cols_unsort, data_sn)
        div = len(meta_isf_rows)  # =len(meta_isf_cols)
        ϵ_st_simple_dict = {}  # 简化版支持用meta_isf_rows的元素拿ϵ_st
        for single_row_st in meta_isf_rows:
            σ_st, μ_st, β_st = single_row_st if len(single_row_st) == 3 else (*single_row_st, None)
            flag, ϵ_st_dict = load_ϵ(self.s_t, σ_st, μ_st, meta_ν_st, β_st, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_ϵ fail by self.s_t={}, σ_st={}, μ_st={}, meta_ν_st={}, β_st={} " \
                          "with msg={}".format(self.s_t, σ_st, μ_st, meta_ν_st, β_st, ϵ_st_dict)
                logger.error(err_msg)
                return False, err_msg
            ϵ_st_simple_dict[single_row_st] = ϵ_st_dict
        ϵ_simple_dict = {}  # 简化版支持用meta_isf_cols的元素拿ϵ
        for single_col in meta_isf_cols:
            ν, β = single_col if isinstance(single_col, tuple) else (single_col, None)
            ϵ_dict = data_σ_μ.get_ϵ_dict(ν, β)
            ϵ_simple_dict[single_col] = ϵ_dict

        # 计算7个对称ISF
        for ϵ_mode in ["ϵ_1", "ϵ_4", "ϵ_14", "ϵ_5", "ϵ_15", "ϵ_6", "ϵ_16"]:
            if ϵ_mode == "ϵ_1":
                # 2.1.1，根据ϵ1算\\[μ][σ][ν']\\
                ϵ_tuple = (self.s_n, data_σ_μ.μ, data_σ_μ.σ, meta_ν_st)
            elif ϵ_mode == "ϵ_4":
                # 2.1.2，根据ϵ4算\\[σ~][μ~][ν‘]\\
                ϵ_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.σ), data_sn.get_dagger(data_σ_μ.μ), meta_ν_st)
            elif ϵ_mode == "ϵ_14":
                # 2.1.3，根据ϵ14算\\[μ~][σ~][ν’]\\
                ϵ_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.μ), data_sn.get_dagger(data_σ_μ.σ), meta_ν_st)
            elif ϵ_mode == "ϵ_5":
                # 2.1.4，根据ϵ5算\\[σ~][μ][ν‘~]\\
                ϵ_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.σ), data_σ_μ.μ, data_st.get_dagger(meta_ν_st))
            elif ϵ_mode == "ϵ_15":
                # 2.1.5，根据ϵ15算\\[μ][σ~][ν’~]\\
                ϵ_tuple = (self.s_n, data_σ_μ.μ, data_sn.get_dagger(data_σ_μ.σ), data_st.get_dagger(meta_ν_st))
            elif ϵ_mode == "ϵ_6":
                # 2.1.6，根据ϵ6算\\[σ][μ~][ν‘~]\\
                ϵ_tuple = (self.s_n, data_σ_μ.σ, data_sn.get_dagger(data_σ_μ.μ), data_st.get_dagger(meta_ν_st))
            else:
                # 2.1.7，根据ϵ16算\\[μ~][σ][ν’~]\\
                ϵ_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.μ), data_σ_μ.σ, data_st.get_dagger(meta_ν_st))

            _, is_calc_ed = is_isf_exist(*ϵ_tuple)  # 为了图方便，这里直接用统一的is_isf_exist代替详细检查对称性了
            if is_calc_ed is False:
                isf_start_time = time.time()

                if ϵ_mode == "ϵ_1":
                    # 2.1.1，[μ'][σ']β'//[ν]β
                    symmetry_isf_rows_unsort = [(i[1], i[0]) if len(i) == 2 else (i[1], i[0], i[2]) for i in
                                                meta_isf_rows]
                elif ϵ_mode == "ϵ_4":
                    # 2.1.2，[σ'~][μ'~]β'//[ν]β
                    symmetry_isf_rows_unsort = [(data_st.get_dagger(i[0]), data_st.get_dagger(i[1])) if len(i) == 2 else
                                                (data_st.get_dagger(i[0]), data_st.get_dagger(i[1]), i[2])
                                                for i in meta_isf_rows]
                elif ϵ_mode == "ϵ_14":
                    # 2.1.3，[μ'~][σ'~]β'//[ν]β
                    symmetry_isf_rows_unsort = [(data_st.get_dagger(i[1]), data_st.get_dagger(i[0])) if len(i) == 2 else
                                                (data_st.get_dagger(i[1]), data_st.get_dagger(i[0]), i[2])
                                                for i in meta_isf_rows]
                elif ϵ_mode == "ϵ_5":
                    # 2.1.4，[σ'~][μ']β'//[ν~]β
                    symmetry_isf_rows_unsort = [(data_st.get_dagger(i[0]), i[1]) if len(i) == 2 else
                                                (data_st.get_dagger(i[0]), i[1], i[2])
                                                for i in meta_isf_rows]
                elif ϵ_mode == "ϵ_15":
                    # 2.1.5，[μ'][σ'~]β'//[ν~]β
                    symmetry_isf_rows_unsort = [(i[1], i[0]) if len(i) == 2 else (i[1], i[0], i[2]) for i in
                                                meta_isf_rows]
                elif ϵ_mode == "ϵ_6":
                    # 2.1.6，[σ'][μ'~]β'//[ν~]β
                    symmetry_isf_rows_unsort = [(i[1], i[0]) if len(i) == 2 else (i[1], i[0], i[2]) for i in
                                                meta_isf_rows]
                else:
                    # 2.1.7，[μ'~][σ']β'//[ν~]β
                    symmetry_isf_rows_unsort = [(i[1], i[0]) if len(i) == 2 else (i[1], i[0], i[2]) for i in
                                                meta_isf_rows]

                symmetry_isf_rows, meta_2_sym_relationship_rows, _ = self._sort_isf_rows(symmetry_isf_rows_unsort,
                                                                                         data_st)
                symmetry_isf_cols = meta_isf_cols if ϵ_mode in ["ϵ_1", "ϵ_4", "ϵ_14"] else dagger_isf_cols_unsort
                meta_2_sym_relationship_cols = list(range(div)) if ϵ_mode in ["ϵ_1", "ϵ_4", "ϵ_14"] else \
                    dagger_meta_2_sym_relationship_cols
                symmetry_isf_square_matrix = \
                    self._calc_symmetry_isf_square_matrix(meta_isf_square_matrix, meta_isf_rows, meta_isf_cols, div,
                                                          ϵ_st_simple_dict, ϵ_simple_dict, ϵ_mode,
                                                          meta_2_sym_relationship_rows, meta_2_sym_relationship_cols)
                symmetry_isf_square = {"rows": symmetry_isf_rows,
                                       "cols": symmetry_isf_cols,
                                       "isf": symmetry_isf_square_matrix}
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
                # 实有β的才会是real list，没有的会是[]
                β_st_list = [i[2] for i in symmetry_isf_rows_unsort if len(i) == 3 and i[0] == σ_st and i[1] == μ_st]
                if β_st_list:
                    for β_st in sorted(β_st_list):  # 排序铜牌的β
                        single_tuple = (σ_st, μ_st, β_st)
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
    def _sort_isf_cols(dagger_isf_cols_unsort, data_sn):
        """将未按照Yamanouchi排序的dagger_isf_cols_unsort按照Yamanouchi排序
        技巧是利用yd_index达到排序yd的目的"""
        sort_isf_cols = []
        sym_2_meta_relationship = []  # 形如[1, 3, 2]表示已排序的cols是从meta中哪里来的
        all_yd_list = data_sn.yd_list
        ν_index_list = [all_yd_list.index(i[0]) if isinstance(i, tuple) else all_yd_list.index(i)
                        for i in dagger_isf_cols_unsort]
        ν_index_one_list = list(set(ν_index_list))
        for ν_index_one in sorted(ν_index_one_list):  # 排序金牌的ν
            ν = all_yd_list[ν_index_one]
            # 实有β的才会是real list，没有的会是[]
            β_list = [i[1] for i in dagger_isf_cols_unsort if isinstance(i, tuple) and i[0] == ν]
            if β_list:
                for β in sorted(β_list):  # 排序银牌的β
                    single_tuple = (ν, β)
                    sort_isf_cols.append(single_tuple)
                    sym_2_meta_relationship.append(dagger_isf_cols_unsort.index(single_tuple))
            else:
                sort_isf_cols.append(ν)
                sym_2_meta_relationship.append(dagger_isf_cols_unsort.index(ν))
        # 形如[1, 3, 2]表示meta会到已排序的cols哪里去
        meta_2_sym_relationship = [sym_2_meta_relationship.index(i) for i in range(len(sym_2_meta_relationship))]
        return sort_isf_cols, meta_2_sym_relationship, sym_2_meta_relationship

    @staticmethod
    def _calc_symmetry_isf_square_matrix(meta_isf_square_matrix, meta_isf_rows, meta_isf_cols, div,
                                         ϵ_st_simple_dict, ϵ_simple_dict, ϵ_key,
                                         meta_2_sym_relationship_rows, meta_2_sym_relationship_cols):
        symmetry_isf_square_matrix = sp.zeros(div)
        for (row_i, row_st), (col_i, col) in product(enumerate(meta_isf_rows), enumerate(meta_isf_cols)):
            sym_row, sym_col = meta_2_sym_relationship_rows[row_i], meta_2_sym_relationship_cols[col_i]
            ϵ = ϵ_simple_dict[col][ϵ_key] * ϵ_st_simple_dict[row_st][ϵ_key]
            symmetry_isf_square_matrix[sym_row, sym_col] = ϵ * meta_isf_square_matrix[row_i, col_i]
        return symmetry_isf_square_matrix


class CGCHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(CGCHelper, self).__init__()

    def calc_meta_cgc_and_ϵ_include_save(self, meta_isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
        """
        凭借元ISF计算元CGC，并且，计算ϵ
        存储元CGC以及ϵ_dict

        1.2.1，判断，并注册元
        1.2.2，计算元CGC
        1.2.3，计算7个ϵ
        """
        # 准备数据
        # 拆解isf_square_dict
        σ_μ_β_all_st_tuple_list = meta_isf_square_dict["rows"]
        ν_β_list = meta_isf_square_dict["cols"]
        isf_square_matrix = meta_isf_square_dict["isf"]
        # 得到h_ν_st
        h_ν_st = data_st.get_yt_num(ν_st)
        #
        is_σ_self_conjugate = (data_σ_μ.σ == data_sn.get_dagger(data_σ_μ.σ))
        is_μ_self_conjugate = (data_σ_μ.μ == data_sn.get_dagger(data_σ_μ.μ))
        is_ν_st_self_conjugate = (ν_st == data_sn.get_dagger(ν_st))
        # 同表自对称条件
        self_symmetry_part = (is_σ_self_conjugate or is_μ_self_conjugate) and is_ν_st_self_conjugate

        # 下面的循环，由互不相干的三层构成：行、列、m。更改它们的顺序不影响计算结果。当前使用，先m再列最后行的顺序，比较优
        for m_ν_st in range(1, h_ν_st + 1):
            # 在深层循环前，准备数据，空间换时间
            cgc_st_square_dict_list_by_row = []
            offset_σ_list = []
            offset_μ_list = []
            for σ_μ_β_all_st in σ_μ_β_all_st_tuple_list:
                σ_st, μ_st, β_st = σ_μ_β_all_st if len(σ_μ_β_all_st) == 3 else (*σ_μ_β_all_st, None)
                flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
                                                    is_flag_true_if_not_s_n=False)
                if not flag:
                    err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, " \
                              "m_ν_st={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
                                                         cgc_st_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
                cgc_st_square_dict_list_by_row.append(cgc_st_square_dict)
                offset_σ = data_sn.get_offset(σ_st, data_σ_μ.σ)
                offset_μ = data_sn.get_offset(μ_st, data_σ_μ.μ)
                offset_σ_list.append(offset_σ)
                offset_μ_list.append(offset_μ)

            for col, ν_β in enumerate(ν_β_list):  # 按照列循环
                # 1.2，根据元ISF计算元CGC【[σ][μ][ν]β】
                # 拆解列信息
                isf_square_vector = isf_square_matrix.col(col)
                ν, β = ν_β if isinstance(ν_β, tuple) else (ν_β, None)
                # 每个ν有自己的offset
                offset_ν = data_sn.get_offset(ν_st, ν)
                m_ν = offset_ν + m_ν_st
                '''
                是元ISF，但不是元CGC的情况
                a，同ISF表，自对称，ev<0的不是元CGC
                b，元ISF，但部分列，是对称的，不是元CGC
                '''
                # a，同ISF表，自对称
                if self_symmetry_part and data_sn.get_eigenvalue(ν) < 0:
                    continue
                # b，元ISF，但部分列，是对称的，不是元CGC
                _, is_calc_ed = is_cgc_exist(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν)
                if is_calc_ed is True:
                    continue
                # 1.2.1，判断，并注册元信息
                data_σ_μ.register_meta_cgc_ν_list(ν, β)

                # 1.2.2，已经得到σ,μ,ν,β,m，开始计算元CGC
                cgc_start_time = time.time()
                meta_cgc_square_dict = {}
                for σ_μ_β_all_st, isf_square, cgc_st_square_dict, offset_σ, offset_μ in \
                        zip(σ_μ_β_all_st_tuple_list, isf_square_vector, cgc_st_square_dict_list_by_row,
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
                              "ν={}, β={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "meet cgc_square_n={} not eq 1".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν,
                                                                     meta_cgc_square_dict, cgc_square_n)
                    logger.error(err_msg)
                    return False, err_msg

                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, meta_cgc_square_dict, cgc_speed_time)
                if not flag:
                    err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, β={}, m_ν={}, meta_cgc_square_dict={} with " \
                              "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, meta_cgc_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

                # 1.2.3，计算7个ϵ
                h_ν = data_sn.get_yt_num(ν)
                is_ν_self_conjugate = (ν == data_sn.get_dagger(ν))
                # if m_ν == 1 and (β is None or β == 1):
                if m_ν == 1:
                    ϵ_dict = self.calc_ϵ_by_m_1(meta_cgc_square_dict, ν, data_sn, data_σ_μ,
                                                is_σ_self_conjugate, is_μ_self_conjugate, is_ν_self_conjugate)
                    data_σ_μ.register_ϵ_dict(ν, β, ϵ_dict)
                    if len(ϵ_dict) == 7:  # ϵ满状态有7个
                        flag, msg = save_ϵ(self.s_n, data_σ_μ.μ, data_σ_μ.σ, ν, β, ϵ_dict)
                        if not flag:
                            err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.μ={}, data_σ_μ.σ={}, " \
                                      "ν={}, β={}, ϵ_dict={} with " \
                                      "msg={}".format(self.s_n, data_σ_μ.μ, data_σ_μ.σ, ν, β, ϵ_dict, msg)
                            logger.error(err_msg)
                            return False, err_msg

                # if m_ν == h_ν and (β is None or β == 1):
                if m_ν == h_ν:
                    ϵ_dict = data_σ_μ.get_ϵ_dict(ν, β)
                    if len(ϵ_dict) == 7:  # ϵ满状态有7个
                        continue
                    adding_ϵ_dict = self.calc_ϵ_by_m_h(meta_cgc_square_dict, ν, data_sn, data_σ_μ,
                                                       is_σ_self_conjugate, is_μ_self_conjugate, is_ν_self_conjugate)
                    adding_ϵ_dict.update(ϵ_dict)
                    data_σ_μ.register_ϵ_dict(ν, β, adding_ϵ_dict)
                    if len(ϵ_dict) == 7:  # ϵ满状态有7个
                        flag, msg = save_ϵ(self.s_n, data_σ_μ.μ, data_σ_μ.σ, ν, β, adding_ϵ_dict)
                        if not flag:
                            err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.μ={}, data_σ_μ.σ={}, " \
                                      "ν={}, β={}, adding_ϵ_dict={} with " \
                                      "msg={}".format(self.s_n, data_σ_μ.μ, data_σ_μ.σ, ν, β, adding_ϵ_dict, msg)
                            logger.error(err_msg)
                            return False, err_msg

        return True, None

    def calc_symmetry_cgc_dict_include_save(self, ν, β, data_sn, data_st, data_σ_μ):
        """
        利用元CGC和ϵ计算全部对称CGC

        3.1.1，load元CGC
        3.2.1，根据ϵ1算【[μ][σ][ν]β】
        3.2.2，根据ϵ4算【[σ~][μ~][ν]β】
        3.2.3，根据ϵ14算【[μ~][σ~][ν]β】
        3.2.4，根据ϵ5算【[σ~][μ][ν~]β】
        3.2.5，根据ϵ15算【[μ][σ~][ν~]β】
        3.2.6，根据ϵ6算【[σ][μ~][ν~]β】
        3.2.7，根据ϵ16算【[μ~][σ][ν~]β】
        """
        h_ν = data_sn.get_yt_num(ν)
        ϵ_dict = data_σ_μ.get_ϵ_dict(ν, β)
        for m in range(1, h_ν + 1):
            m_dagger = h_ν + 1 - m
            # 3.1.1，load元CGC
            flag, meta_cgc_square_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m,
                                                  is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, ν={}, β={}, m={}, " \
                          "with msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m, meta_cgc_square_dict)
                logger.error(err_msg)
                return False, err_msg

            # 一起算所有CGC太占内存了，所以还是分开算，反正最后的优化方向是这里不算，使用时再现算
            # 3.2.1，根据ϵ1算【[μ][σ][ν]β m】
            cgc_start_time = time.time()
            ϵ1_tuple = (self.s_n, data_σ_μ.μ, data_σ_μ.σ, ν, β, m)
            _, is_calc_ed = is_cgc_exist(*ϵ1_tuple)
            if is_calc_ed is False:  # 能进行到这里的ϵ_x必须非None
                ϵ_1 = ϵ_dict["ϵ_1"]
                cgc1_square_dict = self.calc_symmetry_cgc_by_ϵ_1(ϵ_1, meta_cgc_square_dict)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ1_tuple, cgc1_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ1_tuple={}, cgc1_square_dict={} with " \
                              "msg={}".format(ϵ1_tuple, cgc1_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.2，根据ϵ4算【[σ~][μ~][ν]β m】
            cgc_start_time = time.time()
            ϵ4_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.σ), data_sn.get_dagger(data_σ_μ.μ), ν, β, m)
            _, is_calc_ed = is_cgc_exist(*ϵ4_tuple)
            if is_calc_ed is False:
                ϵ_4 = ϵ_dict["ϵ_4"]
                cgc4_square_dict = self._calc_symmetry_cgc_by_ϵ_4(ϵ_4, meta_cgc_square_dict, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ4_tuple, cgc4_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ4_tuple={}, cgc4_square_dict={} with " \
                              "msg={}".format(ϵ4_tuple, cgc4_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.3，根据ϵ14算【[μ~][σ~][ν]β m】
            cgc_start_time = time.time()
            ϵ14_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.μ), data_sn.get_dagger(data_σ_μ.σ), ν, β, m)
            _, is_calc_ed = is_cgc_exist(*ϵ14_tuple)
            if is_calc_ed is False:
                ϵ_14 = ϵ_dict["ϵ_14"]
                cgc14_square_dict = self._calc_symmetry_cgc_by_ϵ_14(ϵ_14, meta_cgc_square_dict, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ14_tuple, cgc14_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ14_tuple={}, cgc14_square_dict={} with " \
                              "msg={}".format(ϵ14_tuple, cgc14_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.4，根据ϵ5算【[σ~][μ][ν~]β m~】
            cgc_start_time = time.time()
            ϵ5_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.σ), data_σ_μ.μ, data_sn.get_dagger(ν), β, m_dagger)
            _, is_calc_ed = is_cgc_exist(*ϵ5_tuple)
            if is_calc_ed is False:
                ϵ_5 = ϵ_dict.get("ϵ_5")
                cgc5_square_dict = self._calc_symmetry_cgc_by_ϵ_5(ϵ_5, meta_cgc_square_dict, ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ5_tuple, cgc5_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ5_tuple={}, cgc5_square_dict={} with " \
                              "msg={}".format(ϵ5_tuple, cgc5_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.5，根据ϵ15算【[μ][σ~][ν~]β m~】
            cgc_start_time = time.time()
            ϵ15_tuple = (self.s_n, data_σ_μ.μ, data_sn.get_dagger(data_σ_μ.σ), data_sn.get_dagger(ν), β, m_dagger)
            _, is_calc_ed = is_cgc_exist(*ϵ15_tuple)
            if is_calc_ed is False:
                ϵ_15 = ϵ_dict.get("ϵ_15")
                cgc15_square_dict = self._calc_symmetry_cgc_by_ϵ_15(ϵ_15, meta_cgc_square_dict, ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ15_tuple, cgc15_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ15_tuple={}, cgc15_square_dict={} with " \
                              "msg={}".format(ϵ15_tuple, cgc15_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.6，根据ϵ6算【[σ][μ~][ν~]β m~】
            cgc_start_time = time.time()
            ϵ6_tuple = (self.s_n, data_σ_μ.σ, data_sn.get_dagger(data_σ_μ.μ), data_sn.get_dagger(ν), β, m_dagger)
            _, is_calc_ed = is_cgc_exist(*ϵ6_tuple)
            if is_calc_ed is False:
                ϵ_6 = ϵ_dict.get("ϵ_6")
                cgc6_square_dict = self._calc_symmetry_cgc_by_ϵ_6(ϵ_6, meta_cgc_square_dict, ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ6_tuple, cgc6_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ6_tuple={}, cgc6_square_dict={} with " \
                              "msg={}".format(ϵ6_tuple, cgc6_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

            # 3.2.7，根据ϵ16算【[μ~][σ][ν~]β m~】
            cgc_start_time = time.time()
            ϵ16_tuple = (self.s_n, data_sn.get_dagger(data_σ_μ.μ), data_σ_μ.σ, data_sn.get_dagger(ν), β, m_dagger)
            _, is_calc_ed = is_cgc_exist(*ϵ16_tuple)
            if is_calc_ed is False:
                ϵ_16 = ϵ_dict.get("ϵ_16")
                cgc16_square_dict = self._calc_symmetry_cgc_by_ϵ_16(ϵ_16, meta_cgc_square_dict, ν, m, data_sn, data_σ_μ)
                cgc_speed_time = int(time.time() - cgc_start_time)
                flag, msg = save_cgc(*(*ϵ16_tuple, cgc16_square_dict, cgc_speed_time))
                if not flag:
                    err_msg = "save_cgc fail by ϵ16_tuple={}, cgc16_square_dict={} with " \
                              "msg={}".format(ϵ16_tuple, cgc16_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

        return True, None

    '''
    在以下计算中，需要注意变量的归属问题。
    使用data_σ_μ得到的变量，说明它的计算中，参数是CGC_σμνβ的
    而使用ϵ_xxx_tuple的，说明它的计算中，参数是对称目标的
    
    另外，用行列式储存CGC，能得到极大的计算速度和简洁性优化；
    但是，更好的优化是根本不计算它们，只留一个索引，需要调取时，根据索引现场计算，来节省计算和储存资源
    '''
    @staticmethod
    def calc_ϵ_by_m_1(cgc_square_dict, ν, data_sn, data_σ_μ,
                      is_σ_self_conjugate, is_μ_self_conjugate, is_ν_self_conjugate):
        """根据元CGC的m_1计算ϵ

        注意：输入的cgc_square_dict的m应该等于1，它需要在外部确认
        ϵ1：公式4-116b
        ϵ4：公式4-122a
        ϵ5：公式4-122b
        ϵ6：公式4-123
        ϵ14：自行推导，使ϵ1作用在经过ϵ4对称后的CGC上
        ϵ15：自行推导，使ϵ1作用在经过ϵ5对称后的CGC上
        ϵ16：自行推导，使ϵ1作用在经过ϵ6对称后的CGC上
        None 值表示形式ϵ，因为实际计算中不需要使用，元CGC早已算完
        """
        ϵ_dict = {}

        # 公共条件
        phase_vector_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_vector_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        # phase_vector_list_ν = data_sn.get_phase_factor_list(ν)
        # Λ_m_ν_1 = phase_vector_list_ν[0]  # =1
        # Λ_h_ν = phase_vector_list_ν[-1]
        # h_σ = data_sn.get_yt_num(data_σ_μ.σ)
        # h_μ = data_sn.get_yt_num(data_σ_μ.μ)
        # TODO 真想把字典换成行列式啊
        m_σ_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple))
        m_μ_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple))

        if data_σ_μ.σ == data_σ_μ.μ:
            ϵ_1 = None
            ϵ_dict["ϵ_1"] = ϵ_1
        else:
            # ϵ1条件【[μ][σ][ν]β】：m_v=1，(m_μ, m_σ)=min，m取（m_σ, m_μ, 1） # 注意它是m_μ在前
            min_m_μ = min(m_μ_set)
            m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[1] == min_m_μ)
            min_m_σ = min(m_σ_second_set)
            ϵ_1_key = (min_m_σ, min_m_μ,)
            ϵ_1 = sp.sign(cgc_square_dict[ϵ_1_key])
            ϵ_dict["ϵ_1"] = ϵ_1

        if is_σ_self_conjugate and is_μ_self_conjugate:
            ϵ_4 = None
            ϵ_dict["ϵ_4"] = ϵ_4
            ϵ_14 = None
            ϵ_dict["ϵ_14"] = ϵ_14
        else:
            # ϵ4条件【[σ~][μ~][ν]β】：m_ν=1，(m_σ, m_μ)=max，m取（m_σ, m_μ, 1）
            max_m_σ = max(m_σ_set)
            m_μ_second_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[0] == max_m_σ)
            max_m_μ = max(m_μ_second_set)
            Λ_max_m_σ = phase_vector_list_σ[max_m_σ - 1]
            Λ_max_m_μ = phase_vector_list_μ[max_m_μ - 1]
            ϵ_4_key = (max_m_σ, max_m_μ,)
            ϵ_4 = sp.sign(cgc_square_dict[ϵ_4_key] * Λ_max_m_σ * Λ_max_m_μ)
            ϵ_dict["ϵ_4"] = ϵ_4

            # ϵ14条件【[μ~][σ~][ν]β】：m_ν=1，(m_μ, m_σ)=max，m取（m_σ, m_μ, 1）
            max_m_μ = max(m_μ_set)
            m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[1] == max_m_μ)
            max_m_σ = max(m_σ_second_set)
            Λ_max_m_μ = phase_vector_list_μ[max_m_μ - 1]
            Λ_max_m_σ = phase_vector_list_σ[max_m_σ - 1]
            ϵ_14_key = (max_m_σ, max_m_μ,)
            ϵ_14 = sp.sign(cgc_square_dict[ϵ_14_key] * Λ_max_m_μ * Λ_max_m_σ)
            ϵ_dict["ϵ_14"] = ϵ_14

        # 对于判断条件，满足它就使得[σ~][μ][ν~]β的实际m_h就是[σ][μ][ν]β的m1了。
        # 注意，此时ϵ5和ϵ15都只是形式上的了，实际计算中不需要使用，因为元CGC早已算完
        if is_σ_self_conjugate and is_ν_self_conjugate:  # 有了它就使得[σ~][μ][ν~]β的实际m1就是[σ][μ][ν]β的m1了
            # 此时【[σ~][μ][ν~]β】=【[σ][μ][ν]β】(但m不一致)
            ϵ_5 = None
            ϵ_dict["ϵ_5"] = ϵ_5

            # 此时【[μ][σ~][ν~]β】=【[μ][σ][ν]β】(但m不一致)
            ϵ_15 = None
            ϵ_dict["ϵ_15"] = ϵ_15

        # 对于判断条件，满足它就使得[σ][μ~][ν~]β的实际m1就是[σ][μ][ν]β的m1了。
        # 注意，此时ϵ6和ϵ16都只是形式上的了，实际计算中不需要使用，因为元CGC早已算完
        if is_μ_self_conjugate and is_ν_self_conjugate:
            # 此时【[σ][μ~][ν~]β】=【[σ][μ][ν]β】(但m不一致)
            ϵ_6 = None
            ϵ_dict["ϵ_6"] = ϵ_6

            # 此时【[μ~][σ][ν~]β】=【[μ][σ][ν]β】(但m不一致)
            ϵ_16 = None
            ϵ_dict["ϵ_16"] = ϵ_16

        return ϵ_dict

    @staticmethod
    def calc_ϵ_by_m_h(cgc_square_dict, ν, data_sn, data_σ_μ,
                      is_σ_self_conjugate, is_μ_self_conjugate, is_ν_self_conjugate):
        """根据元CGC的m_h计算ϵ

        注意：输入的cgc_square_dict的m应该等于h_ν，它需要在外部确认
        ϵ5：公式4-122b
        ϵ6：公式4-123
        ϵ15：自行推导，使ϵ1作用在经过ϵ5对称后的CGC上
        ϵ16：自行推导，使ϵ1作用在经过ϵ6对称后的CGC上
        """
        ϵ_h_dict = {}
        if is_σ_self_conjugate and is_μ_self_conjugate and is_ν_self_conjugate:
            # 这种情况不用算，快速退出
            return ϵ_h_dict
        # 公共条件
        phase_vector_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_vector_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        phase_vector_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_h_ν = phase_vector_list_ν[-1]
        h_σ = data_sn.get_yt_num(data_σ_μ.σ)
        h_μ = data_sn.get_yt_num(data_σ_μ.μ)
        m_σ_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple))
        m_μ_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple))

        if not is_σ_self_conjugate or not is_ν_self_conjugate:
            # ϵ5条件【[σ~][μ][ν~]β】：m_ν=h_ν，(m_σ, m_μ)=min，m取（h_σ+1-m_σ, m_μ, h_ν）
            min_m_σ = min(m_σ_set)
            m_μ_second_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[0] == min_m_σ)
            min_m_μ = min(m_μ_second_set)
            Λ_min_m_σ = phase_vector_list_σ[min_m_σ - 1]  # 物理序要转换成py序
            ϵ_5_key = (h_σ + 1 - min_m_σ, min_m_μ,)
            ϵ_5 = sp.sign(cgc_square_dict[ϵ_5_key] * Λ_min_m_σ * Λ_h_ν)
            ϵ_h_dict["ϵ_5"] = ϵ_5

            # ϵ15条件【[μ][σ~][ν~]β】：m_ν=h_ν，(m_μ, m_σ)=min，m取（h_σ+1-m_σ, m_μ, h_ν）
            min_m_μ = min(m_μ_set)
            m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[1] == min_m_μ)
            min_m_σ = min(m_σ_second_set)
            Λ_min_m_σ = phase_vector_list_σ[min_m_σ - 1]  # 物理序要转换成py序
            Λ_h_ν = phase_vector_list_ν[-1]
            ϵ_15_key = (h_σ + 1 - min_m_σ, min_m_μ,)
            ϵ_15 = sp.sign(cgc_square_dict[ϵ_15_key] * Λ_min_m_σ * Λ_h_ν)
            ϵ_h_dict["ϵ_15"] = ϵ_15

        if not is_μ_self_conjugate or not is_ν_self_conjugate:
            # ϵ6条件【[σ][μ~][ν~]β】：m_ν=h_ν，(m_σ, m_μ)=min，m取（m_σ, h_μ+1-m_μ, h_ν）
            min_m_σ = min(m_σ_set)
            m_μ_second_set = set(i[1] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[0] == min_m_σ)
            min_m_μ = min(m_μ_second_set)
            Λ_min_m_μ = phase_vector_list_μ[min_m_μ - 1]  # 物理序要转换成py序
            ϵ_6_key = (min_m_σ, h_μ + 1 - min_m_μ,)
            ϵ_6 = sp.sign(cgc_square_dict[ϵ_6_key] * Λ_min_m_μ * Λ_h_ν)
            ϵ_h_dict["ϵ_6"] = ϵ_6

            # ϵ16条件【[μ~][σ][ν~]β】：m_ν=h_ν，(m_μ, m_σ)=min，m取（m_σ, h_μ+1-m_μ, h_ν）
            min_m_μ = min(m_μ_set)
            m_σ_second_set = set(i[0] for i in cgc_square_dict.keys() if isinstance(i, tuple) and i[1] == min_m_μ)
            min_m_σ = min(m_σ_second_set)
            Λ_min_m_μ = phase_vector_list_μ[min_m_μ - 1]  # 物理序要转换成py序
            ϵ_16_key = (min_m_σ, h_μ + 1 - min_m_μ,)
            ϵ_16 = sp.sign(cgc_square_dict[ϵ_16_key] * Λ_min_m_μ * Λ_h_ν)
            ϵ_h_dict["ϵ_16"] = ϵ_16

        return ϵ_h_dict

    '''下面的calc_symmetry_cgc_by_ϵ_xxx不怕重复，因为验证过后，它们将不会被计算'''

    @staticmethod
    def calc_symmetry_cgc_by_ϵ_1(ϵ_1, meta_cgc_square_dict):
        # 根据ϵ1算【[μ][σ][ν]β m】
        cgc1_square_dict = {(k[1], k[0]): ϵ_1 * v for k, v in meta_cgc_square_dict.items()}
        return cgc1_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_4(ϵ_4, meta_cgc_square_dict, data_sn, data_σ_μ):
        # 根据ϵ4算【[σ~][μ~][ν]β m】
        cgc4_square_dict = {}
        h_σ_add_1, h_μ_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(data_σ_μ.μ) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc4_square_dict[(h_σ_add_1 - m_σ, h_μ_add_1 - m_μ)] = ϵ_4 * Λ_m_σ * Λ_m_μ * v
        return cgc4_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_14(ϵ_14, meta_cgc_square_dict, data_sn, data_σ_μ):
        # 根据ϵ14算【[μ~][σ~][ν]β m】
        cgc14_square_dict = {}
        h_σ_add_1, h_μ_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(data_σ_μ.μ) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc14_square_dict[(h_μ_add_1 - m_μ, h_σ_add_1 - m_σ)] = ϵ_14 * Λ_m_σ * Λ_m_μ * v
        return cgc14_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_5(ϵ_5, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ5算【[σ~][μ][ν~]β m~】
        cgc5_square_dict = {}
        h_σ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            cgc5_square_dict[(h_σ_add_1 - m_σ, m_μ)] = ϵ_5 * Λ_m_σ * Λ_m_ν * v
        return cgc5_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_15(ϵ_15, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ15算【[μ][σ~][ν~]β m~】
        cgc15_square_dict = {}
        h_σ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.σ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_σ = data_sn.get_phase_factor_list(data_σ_μ.σ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_σ = phase_factor_list_σ[m_σ - 1]
            cgc15_square_dict[(m_μ, h_σ_add_1 - m_σ)] = ϵ_15 * Λ_m_σ * Λ_m_ν * v
        return cgc15_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_6(ϵ_6, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ6算【[σ][μ~][ν~]β m~】
        cgc6_square_dict = {}
        h_μ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.μ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc6_square_dict[(m_σ, h_μ_add_1 - m_μ)] = ϵ_6 * Λ_m_μ * Λ_m_ν * v
        return cgc6_square_dict

    @staticmethod
    def _calc_symmetry_cgc_by_ϵ_16(ϵ_16, meta_cgc_square_dict, ν, m_ν, data_sn, data_σ_μ):
        # 根据ϵ16算【[μ~][σ][ν~]β m~】
        cgc16_square_dict = {}
        h_μ_add_1, h_ν_add_1 = data_sn.get_yt_num(data_σ_μ.μ) + 1, data_sn.get_yt_num(ν) + 1
        phase_factor_list_μ = data_sn.get_phase_factor_list(data_σ_μ.μ)
        phase_factor_list_ν = data_sn.get_phase_factor_list(ν)
        Λ_m_ν = phase_factor_list_ν[m_ν - 1]
        for (m_σ, m_μ), v in meta_cgc_square_dict.items():
            Λ_m_μ = phase_factor_list_μ[m_μ - 1]
            cgc16_square_dict[(h_μ_add_1 - m_μ, m_σ)] = ϵ_16 * Λ_m_μ * Λ_m_ν * v
        return cgc16_square_dict

    # def calc_ϵ_and_cgc_dict_by_isf_include_save(self, isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
    #     """
    #     凭借ISF计算CGC，并且，
    #     根据对称性，尽可能多地计算ϵ和对应的CGC；
    #     ϵ_dict和cgc_dict都要储存
    #     """
    #     # 拆解isf_square_dict，准备数据
    #     σ_μ_β_all_st_tuple_list = isf_square_dict["rows"]
    #     ν_β_list = isf_square_dict["cols"]
    #     isf_square_matrix = isf_square_dict["isf"]
    #
    #     # TODO 所以一种可能的优化思路是，把m_ν_st循环放在最外面。这样，这个ISF的行和列都可以先一步提取做成字典，空间换时间
    #
    #     for i, ν_β in enumerate(ν_β_list):  # 按照列，也就是ν+β开启循环
    #         isf_square_vector = isf_square_matrix.col(i)
    #         ν, β = ν_β if isinstance(ν_β, tuple) else (ν_β, None)
    #         h_ν_st = data_st.yt_num_dict[tuple(ν_st)]
    #         offset_of_m_ν = self._calc_m_with_m_st(ν_st, 0, data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)
    #
    #         for m_ν_st in range(1, h_ν_st + 1):
    #             single_cgc_part_start_time = time.time()
    #             m_ν = offset_of_m_ν + m_ν_st
    #
    #             cgc_square_dict = {}
    #             for σ_μ_β_all_st, isf_square in zip(σ_μ_β_all_st_tuple_list, isf_square_vector):
    #                 # TODO 这个循环里的load_cgc肯定重复了，它和ν_β无关；但是，按m全取的话，内存压力太大了
    #                 σ_st, μ_st, β_st = σ_μ_β_all_st if len(σ_μ_β_all_st) == 3 else (*σ_μ_β_all_st, None)
    #                 flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
    #                                                     is_flag_true_if_not_s_n=False)
    #                 if not flag:
    #                     err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, " \
    #                               "m_ν_st={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
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
    #                           "ν={}, β={}, m_ν={}, cgc_square_dict={} with " \
    #                           "meet cgc_square_n={} not eq 1".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν,
    #                                                                  cgc_square_dict, cgc_square_n)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #             single_cgc_part_speed_time = int(time.time() - single_cgc_part_start_time)
    #             flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_dict,
    #                                  single_cgc_part_speed_time)
    #             if not flag:
    #                 err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
    #                           "ν={}, β={}, m_ν={}, cgc_square_dict={} with " \
    #                           "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_dict, msg)
    #                 logger.error(err_msg)
    #                 return False, err_msg
    #
    #     return True, None
