# -*- coding: utf-8 -*-
"""
this code for creating ISF and CGC
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）本段代码主要依据置换群CGC的第二种递推计算方法
# ISF：第四章第19节 公式：4-192（用Sn-1的本征方程计算） 4-189bc（用Sn的CG系数计算）
# CGC：第四章第11、12节
# 整体思路是，先求Sn的ISF，再利用Sn的ISF以及Sn-1的CG，拼出Sn的CG


import copy
import json
import time
# import numpy as np
import sympy as sp
from functools import lru_cache, partial
from itertools import product, combinations
from conf.cgc_config import default_s_n, isf_0_error_value
from conf.cgc_config import min_s_n_of_isf, min_s_n_of_cgc, min_s_n_of_branching_law
from core.young_diagrams import load_young_diagrams
from core.young_tableaux import load_young_table_num
from core.branching_laws import load_branching_law
from core.yamanouchi_matrix import load_yamanouchi_matrix
from core.cg_series import load_cg_series
from core.eigenvalues import load_eigenvalues
from core.cgc_utils.cgc_db_typing import ISFInfo, CGCInfo
from core.cgc_utils.cgc_local_db import get_isf_file_name, get_isf_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_cgc_file_name, get_cgc_finish_s_n_name
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
    '''
    循环内命名规范：
    |σ>   σ/yds_σ: 表示Si的单个yd/全部yds
    |σ’>  σ_st/yds_σ_st: 表示St的单个yd/全部yds
    |σ’>  bl_yd_of_σ/bl_yds_of_σ: 表示St的单个yd/全部yds中，符合Si是σ的那些σ_st
    |σ’’> bl_yd_of_σ_st/bl_yds_of_σ_st: 表示St-1的单个yd/全部yds中，符合St是σ_st的那些σ_st_s_t
    上面在不发生歧义的情况下可以去掉of
    cgc_tuple = (σ, μ, ν, β, m): 区分CGC的全部下角标
    isf_tuple = ()
    '''
    data_si = None
    data_st = None
    # 按照从小到大的顺序，逐个计算s_i的eigenvalues并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        s_i_isf_speed_time = 0
        logger.debug("Si={}".format(s_i))

        if s_i == 1:
            # 特殊情况，Sn=1时，没有ISF，只有CGC
            σ, μ, ν, β, m = [1], [1], [1], None, 1
            cgc_square_dict = {(1, 1,): 1, "N": 1}
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

        # Sn循环可以得到的数据(如果可以从上一次循环继承，则继承；否则要计算)
        if data_st is not None:
            data_st_st_yt_num_dict = data_st.yt_num_dict
        elif s_i - 1 - 1 >= min_s_n_of_branching_law:
            data_st_st_yt_num_dict = {}
            _, yd_st_st_list = load_young_diagrams(s_i - 1 - 1)
            for yd_st_st in yd_st_st_list:
                _, total_num = load_young_table_num(s_i - 1 - 1, yd_st_st)
                data_st_st_yt_num_dict[tuple(yd_st_st)] = total_num
        else:
            data_st_st_yt_num_dict = None

        data_st = data_si if data_si is not None else DataHelper(s_i - 1)
        data_si = DataHelper(s_i)

        # σ μ
        for σ, μ in product(data_si.yd_list, repeat=2):  # [σ], [μ]双循环
            # TODO σ, μ = μ, σ的情况，看看是真算用来消减误差，还是不算用来节约算力
            # TODO σ, μ = μ, σ的情况，必须用对称性解决，因为在计算带有自由度的ν时，会由于改变系数顺序，造成高斯消元法得到的结果不对称
            '''
            例如，ISF_1 [3, 2] * [3, 1, 1] = [2, 1, 1] 和 ISF_2 [3, 1, 1] * [3, 2] = [2, 1, 1]
            按照对称性要求，它们的结果只可以有写作约定顺序的不同，不可以是两套旋转了的基失组
            
            ISF_1的isf_matrix是
            Matrix([[1/2, sqrt(5)/2, sqrt(3)/2, sqrt(5)/2], 
                    [sqrt(5)/2, 1/2, -sqrt(15)/2, -1/2], 
                    [sqrt(3)/2, -sqrt(15)/2, -1, 0], 
                    [sqrt(5)/2, -1/2, 0, 1]])
            ISF_2的isf_matrix是
            Matrix([[1/2, -sqrt(3)/2, -sqrt(5)/2, -sqrt(5)/2], 
                    [-sqrt(3)/2, -1, -sqrt(15)/2, 0], 
                    [-sqrt(5)/2, -sqrt(15)/2, 1/2, -1/2], 
                    [-sqrt(5)/2, 0, -1/2, 1]]) 
                    
            如果正常解矩阵，并使用高斯消元法和施密特正交归一化手续，得到的结果分别为
            ISF_1（不违反对称性）
            Matrix([
            [-3/10,    0,   3/5, 1/10],
            [-1/30, 5/12, -3/20,  2/5],
            [ -1/6,  1/3,     0, -1/2],
            [  1/2,  1/4,   1/4,    0]])
            ISF_2(违反了对称性)(注意：不是这种计算有误，只是它不符合我们在自由相位选取上约定的一致性和对称性)
            Matrix([
            [-1/14,  -3/7,  1/3,  1/6],
            [-3/14,  1/28, -1/4,  1/2],
            [ 7/10,     0,    0, 3/10],
            [-1/70, 15/28, 5/12, 1/30]])
            
            所以，应该直接使用对称性，得到按照统一约定自由相位的结果
            ISF_2
            Matrix([
            [   0,  1/2, -1/3,   1/6],
            [ 1/4,    0, -1/4,  -1/2],
            [-3/5, 1/10,    0, -3/10],
            [3/20,  2/5, 5/12, -1/30]])
            详见《群表示论的新途径》4-19节，3中最后一段的论述
            '''
            '''
            σ, μ = μ, σ的情况，
            对于ISF，使用公式4-196b
            既，ISF_σσ'μμ'νν'ββ' = ϵ(σμνβ)ϵ(σ'μ'ν'β') * ISF_μμ'σσ'νν'ββ',
            也既，ISF_μμ'σσ'νν'ββ' = ϵ(σμνβ)ϵ(σ'μ'ν'β') * ISF_σσ'μμ'νν'ββ'（这个式子，可做推倒式子使用）
            对于CGC，使用公式4-116a
            既，CGC_{σm_σ,μm_μ,νβm} = ϵ(σμνβ) * CGC_{μm_μ,σm_σ,νβm}
            也既，CGC_{μm_μ,σm_σ,νβm} = ϵ(σμνβ) * CGC_{σm_σ,μm_μ,νβm} （这个式子，可做推倒式子使用）
            其中，ϵ(σμνβ) = sign(CGC_{σm_σ,μm_μ,νβm_1} | (m_μ,m_σ)=min)  (源自绝对相位约定)
            表示ϵ取(m_μ,m_σ)编号为最小时的非零CG系数的符号(注意，μ在前，σ在后)
            '''
            # σ_μ_start_time = time.time()
            logger.debug("σ={}, μ={}".format(σ, μ))

            # σ, μ循环可以得到的数据
            data_σ_μ = ΣMDataHelper(s_i, σ, μ, data_si, data_st)

            for ν_st in data_st.yd_list:  # [ν']循环
                # # ISF
                single_isf_start_time = time.time()
                # 带着None的rows
                row_index_tmp_list = isf_func.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                                   data_σ_μ.cg_series_st_list_dict,
                                                                   ν_st, data_st.yd_list)
                if not row_index_tmp_list:
                    continue  # σ, μ下没有这个ν_st
                logger.debug("ν_st={}".format(ν_st))
                # isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                #                    "cols":[[ν], ([ν], β), ...],  # 有自由度tuple，无自由度list
                #                    "isf": isf_square_matrix}  # np.array/sp.Matrix([len(rows), len(cols)])
                flag, isf_square_dict = isf_func.calc_isf_dict(row_index_tmp_list, ν_st, data_si, data_st, data_σ_μ,
                                                               data_st_st_yt_num_dict)
                if not flag:
                    err_msg = "calc isf_square_dict meet error by row_index_tmp_list={}, ν_st={}, data_si={}, " \
                              "data_st={}, data_σ_μ={}, data_st_st_yt_num_dict={} " \
                              "with msg={}".format(row_index_tmp_list, ν_st, data_si, data_st, data_σ_μ,
                                                   data_st_st_yt_num_dict, isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg

                if _debug_condition(data_σ_μ, ν_st):
                    logger.warning("@@@@ isf_square_dict={} with σ={}, μ={}, ν_st={}".format(
                        isf_square_dict, data_σ_μ.σ, data_σ_μ.μ, ν_st))

                single_isf_speed_time = int(time.time() - single_isf_start_time)
                s_i_isf_speed_time += single_isf_speed_time
                flag, msg = save_isf(s_i, σ, μ, ν_st, isf_square_dict, single_isf_speed_time)
                if not flag:
                    err_msg = "save CG series meet error with s_i={}, σ={}, μ={}, ν_st={}, isf_square_dict={}, " \
                              "msg={}".format(s_i, σ, μ, ν_st, isf_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

                # # CGC
                # 这里只能部分地拼凑CGC，因为完整的一个CGC需要σ * μ下所有的ν'中，同一个ν的ISF
                flag, msg = cgc_func.calc_cgc_dict_part_and_save_by_isf(isf_square_dict, ν_st,
                                                                        data_si, data_st, data_σ_μ)
                if not flag:
                    err_msg = "calc_cgc_dict_part_and_save fail by ν_st={}, data_si={}, data_st={}, " \
                              "data_σ_μ={} isf_square_dict={} with msg={}".format(ν_st, data_si, data_st, data_σ_μ,
                                                                                  isf_square_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

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
                       "cols": [[ν], ([ν], β), ...],                   # 有自由度tuple，无自由度list
                       "isf": isf_square_matrix}                       # np.array([len(rows), len(cols)], dtype=float)

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
     "isf": np.array([[5/12, 1/2, 1/12, 0],
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


def save_cgc(s_n: int, σ: list, μ: list, ν: list, β: (int, None), m: int, cgc_square_dict: dict, speed_time: int):
    """
    这个db用来存CGC

    与其他保存逻辑不同，cgc是一部分一部分算出来的，所以不仅需要保存，还需要更新

    TODO tmp
    <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_β_m.pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_β_m,
    "data": cgc_square_dict,
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型； (Sn)
    beta对应[ν]的多重性;
    m是[ν]的yamanouchi序;
    cgc_square_dict数值是cgc的平方，符号是平方前cgc的符号;
    speed_time表示计算用时（秒）

    例如，
    <CG>/cgc_info/S4/[2,2]_[3,1]/[2,1,1]_1_m2.pkl
    {(2, 1): 0.4999999999999999, (1, 3): 0.24999999999999994, (2, 2): 0.24999999999999994, 'N': 0.9999999999999998}
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
    if not isinstance(cgc_square_dict, dict):
        err_msg = "cgc_square_dict={} with type={} must be dict".format(cgc_square_dict, type(cgc_square_dict))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = CGCInfo(s_n)
    _, file_name = get_cgc_file_name(s_n, σ, μ, ν, β, m)

    flag, old_data = db_info.query_by_file_name(file_name)
    if not flag:
        err_msg = "save_cgc meet error with file_name={} msg={}".format(file_name, old_data)
        logger.error(err_msg)
        return flag, err_msg

    if old_data is False:  # 没有旧的，说明是新建
        table = {"file_name": file_name,
                 "data": cgc_square_dict,
                 "flags": {"speed_time": speed_time}}
        flag, msg = db_info.insert(table)
        if not flag:
            return flag, msg
        flag, msg = db_info.insert_txt(table)
        if not flag:
            return flag, msg

    else:
        table = {"file_name": file_name,
                 "data": cgc_square_dict,
                 "flags": {"speed_time": speed_time + old_data.get("flags").get("speed_time")}}
        flag, msg = db_info.update_by_file_name(file_name, partial_table=table)
        if not flag:
            return flag, msg
        flag, msg = db_info.update_txt_by_file_name(file_name, partial_table=table, point_key="data")
        if not flag:
            return flag, msg

    return True, None


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


def load_isf(s_n: int, σ: list, μ: list, ν_st: list, is_flag_true_if_not_s_n=True,
             output_mode="", ex_params=None):
    """支持两种load，分别是：

    1，只输入必要参数。返回包含isf矩阵及行、列index的字典；
    2.1，output_mode='single_row'。返回包含单行isf及列index
    2.2，output_mode='single_col'。返回包含单列isf及行index
    2.3，output_mode='single_isf'。返回单独一个isf
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
    # if hasattr(ex_params, "__iter__"):
    #     ex_params_list = list(ex_params)
    #     ex_params_list.remove(None)
    #     ex_params = tuple(ex_params_list)
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
                ex_row = ex_params[: -1] if None in ex_params else ex_params
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
                ex_col = ex_params[0] if None in ex_params else ex_params
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
                ex_row = ex_params[0][: -1] if None in ex_params[0] else ex_params[0]
                ex_col = ex_params[1][0] if None in ex_params[1] else ex_params[1]
                if ex_row in isf_square_dict.get("rows", []) and ex_params[1] in isf_square_dict.get("cols", []):
                    row_index = isf_square_dict.get("rows").index(ex_row)
                    cols_index = isf_square_dict.get("cols").index(ex_col)
                    rst = isf_square_dict["isf"][row_index, cols_index]
                    return True, rst  # bingo(3/4)！
                else:
                    err_msg = "ex_row={} should in rows={} and ex_col={} should in cols={} " \
                              "with isf_square_dict={} but not, pls check".format(
                        ex_row, isf_square_dict.get("rows"),
                        ex_col, isf_square_dict.get("cols"), isf_square_dict)
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
            flag, in_matrix_σ = load_yamanouchi_matrix(self.s_n, self.σ, in_key, mode="in")
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
                flag, in_matrix_μ = load_yamanouchi_matrix(self.s_n, self.μ, in_key, mode="in")
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


class DataHelper(object):
    """这里封装供模块内部，Si循环应该取得的数据"""

    def __init__(self, s_k):
        if not isinstance(s_k, int) or s_k < min(min_s_n_of_isf, min_s_n_of_cgc):
            raise Exception("s_k={} must be int and >= 1".format(s_k))
        self.s_n = s_k  # 既可以实例化Sn，也可以实例化St

        self.yd_list = None
        self.bl_yd_list_dict = None
        self.yt_num_dict = None
        self.eigenvalue_list = None
        self._init_data()

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
        yt_num_dict = {}
        for yd in yd_list:
            flag, total_num = load_young_table_num(self.s_n, yd, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get yd total_num meet error with s_i={}, yd={}, msg={}".format(self.s_n, yd, total_num)
                logger.error(err_msg)
                raise Exception(err_msg)
            yt_num_dict[tuple(yd)] = total_num
        self.yt_num_dict = yt_num_dict

        # 4, eigenvalue，按照yd序
        flag, eigenvalue_list = load_eigenvalues(self.s_n, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get eigenvalue_list meet error with s_n={}, msg={}".format(self.s_n, eigenvalue_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.eigenvalue_list = eigenvalue_list


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

    @staticmethod
    def _calc_m_with_m_st(yd_st, m_st, bl_yd_list, yt_st_num_dict):
        """通过m'计算m，这回，我们可以使用分支律以及yt_num了; 也可以通过m_st=0，计算偏移量"""
        # 是core.young_tableaux.py: quickly_calc_young_table_in_decreasing_page_order的特定情况简化
        bl_yd_index = bl_yd_list.index(yd_st)  # 必须成立
        # 对yd_st的前所有分支的yt_st_num以及本项m_st求和
        m = sum([yt_st_num_dict[tuple(bl_yd)] for bl_yd in bl_yd_list[:bl_yd_index]], m_st)
        return m


class ISFHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(ISFHelper, self).__init__()

    def calc_isf_dict(self, row_index_tmp_list, ν_st, data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict):
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

        flag, isf_matrix = self._calc_isf_matrix(row_index_tmp_list, ν_st, data_st.yt_num_dict, data_σ_μ)
        # TODO 到这里isf_matrix是可以通过符号计算，得到无误差的矩阵的。但后面，求本征值的时候，理论上不一定还能保持符号计算了（一元五次方程无公式解）
        if not flag:
            err_msg = "get isf_matrix by row_index_tmp_list={}, ν_st={}, data_st.yt_num_dict={}, data_σ_μ={} fail " \
                      "with msg={}".format(row_index_tmp_list, ν_st, data_st.yt_num_dict, data_σ_μ, isf_matrix)
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
            eigen_tuple_list = isf_matrix.eigenvects()
        except Exception as e:
            logger.error(Exception(e))
            err_msg = "calc eigh meet fail with isf_matrix={} s_n={}".format(isf_matrix, self.s_n)
            logger.error(err_msg)
            return False, err_msg

        λ_ν_st = data_st.eigenvalue_list[data_st.yd_list.index(ν_st)]

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
                                                data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict)
            if not flag:
                err_msg = "calc _calc_isf_β_and_phase_list meet error " \
                          "with soe_vectors={}, ν={}, ν_st={}, row_index_tmp_list={}, β_max={} " \
                          "data_sn={}, data_st={}, data_σ_μ={}, data_st_st_yt_num_dict={}, " \
                          "msg={}".format(soe_vectors, ν, ν_st, row_index_tmp_list, β_max, data_sn, data_st, data_σ_μ,
                                          data_st_st_yt_num_dict, β_tmp_list)
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
                                   data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict):
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
            m_ν = self._calc_m_with_m_st(ν_st, m_ν_st, data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)
            ν_st_st = data_st.bl_yd_list_dict[tuple(ν_st)][0]
            m_ν_st_fbl = self._calc_m_with_m_st(ν_st_st, m_ν_st, data_st.bl_yd_list_dict[tuple(ν_st_fbl)],
                                                data_st_st_yt_num_dict)
            m_ν_by_m_ν_st_fbl = self._calc_m_with_m_st(ν_st_fbl, m_ν_st_fbl, data_sn.bl_yd_list_dict[tuple(ν)],
                                                       data_st.yt_num_dict)

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
                        row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ)
                    if not flag:
                        err_msg = "calc _calc_isf_fbl_another_way meet error with ν={}, m_ν={}, " \
                                  "m_ν_by_m_ν_st_fbl={}, row={}, ν_st_fbl={}, row_index_tmp_list={}, " \
                                  "soe_vector={}, ν_st={}, m_ν_st={}, data_sn={}, data_st={}, data_σ_μ={}, " \
                                  "msg={}".format(ν, m_ν, m_ν_by_m_ν_st_fbl, row,
                                                  ν_st_fbl, row_index_tmp_list, soe_vector, ν_st, m_ν_st,
                                                  data_sn, data_st, data_σ_μ, isf_fbl_another)
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
                        row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ)
                    if not flag:
                        err_msg = "calc _calc_isf_fbl_another_way meet error with ν={}, m_ν={}, " \
                                  "m_ν_by_m_ν_st_fbl={}, first_no_0_isf_fbl_row={}, ν_st_fbl={}, " \
                                  "row_index_tmp_list={}, soe_vector={}, ν_st={}, m_ν_st={}, " \
                                  "data_sn={}, data_st={}, data_σ_μ={}, " \
                                  "msg={}".format(ν, m_ν, m_ν_by_m_ν_st_fbl, first_no_0_isf_fbl_row,
                                                  ν_st_fbl, row_index_tmp_list, soe_vector, ν_st, m_ν_st,
                                                  data_sn, data_st, data_σ_μ, isf_fbl_another)
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
                                  row_index_tmp_list, soe_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ):
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
        cgc_square_fbl_n = cgc_square_st_fbl_dict.pop("N")

        sum_3_loop = 0
        for (σ_st, μ_st, β_st), soe_vector_element in zip(row_index_tmp_list, soe_vector):
            flag, cgc_square_st_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_st_dict with self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, m_ν_st={} " \
                          "meet error with msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, cgc_square_st_dict)
                logger.error(err_msg)
                return False, err_msg
            cgc_square_dict_n = cgc_square_st_dict.pop("N")
            sum_3_loop_part = 0
            for (m_σ_st_fbl, m_μ_st_fbl), cgc_square_st_fbl_element \
                    in cgc_square_st_fbl_dict.items():
                m_σ_fbl = self._calc_m_with_m_st(σ_st_fbl, m_σ_st_fbl, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
                                                 data_st.yt_num_dict)
                m_μ_fbl = self._calc_m_with_m_st(μ_st_fbl, m_μ_st_fbl, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
                                                 data_st.yt_num_dict)
                for (m_σ_st, m_μ_st), cgc_square_st_element in cgc_square_st_dict.items():
                    m_σ = self._calc_m_with_m_st(σ_st, m_σ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
                                                 data_st.yt_num_dict)
                    m_μ = self._calc_m_with_m_st(μ_st, m_μ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
                                                 data_st.yt_num_dict)
                    in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_fbl - 1, m_σ - 1]
                    in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_fbl - 1, m_μ - 1]
                    sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
                                       * sp.sign(cgc_square_st_element) * sp.sign(cgc_square_st_fbl_element) \
                                       * sp.sqrt(abs(cgc_square_st_element * cgc_square_st_fbl_element
                                                     / (cgc_square_dict_n * cgc_square_fbl_n)))
            sum_3_loop += sum_3_loop_part * soe_vector_element

        flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in")  # (Sn-1, Sn)的对换
        if not flag:
            err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
                      "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
            logger.error(err_msg)
            return False, err_msg
        in_matrix_ν_m_m_element = in_matrix_ν[m_ν - 1, m_ν_by_m_ν_st_fbl - 1]
        fbl_isf = sum_3_loop / in_matrix_ν_m_m_element
        return True, fbl_isf

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
    def calc_row_indexes_tmp(bl_yds_of_σ, bl_yds_of_μ, cg_series_st_list_dict, ν_st, yd_st_list):
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, β'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        row_index_tmp_list = []  # [(bl_σ, bl_μ, β'), ...]
        for bl_σ, bl_μ in product(bl_yds_of_σ, bl_yds_of_μ):
            cg_series_st_list = cg_series_st_list_dict[(tuple(bl_σ), tuple(bl_μ))]  # 必须有
            single_ν_cg_series = cg_series_st_list[yd_st_list.index(ν_st)]
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

    def _cgc_st_2_cgc_m_dict(self, cgc_st_tuple, yt_st_num_dict, data_σ_μ):
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
        cgc_st_square_n = cgc_st_square_dict.pop("N")
        for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
            m_σ = self._calc_m_with_m_st(σ_st, m_σ_st, data_σ_μ.bl_yds_of_σ, yt_st_num_dict)
            m_μ = self._calc_m_with_m_st(μ_st, m_μ_st, data_σ_μ.bl_yds_of_μ, yt_st_num_dict)
            rst_dict[(m_σ, m_μ,)] = cgc_st_square
        cgc_st_square_dict["N"] = cgc_st_square_n
        rst_dict["N"] = cgc_st_square_n

        return True, rst_dict
    
    def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right, yt_st_num_dict, data_σ_μ):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        _, factor_cgc_left_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_left, yt_st_num_dict, data_σ_μ)
        left_n = factor_cgc_left_dict.pop("N")
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = sp.sign(factor_cgc_left) * sp.sqrt(abs(factor_cgc_left))
            left_tmp_dict[left_key] = left_tmp
        factor_cgc_left_dict["N"] = left_n  # 还原字典，避免deepcopy
        # 右
        if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
            right_tmp_dict = left_tmp_dict
            left_right_sqrt_n = left_n
        else:
            _, factor_cgc_right_dict = self._cgc_st_2_cgc_m_dict(cgc_st_tuple_right, yt_st_num_dict, data_σ_μ)
            right_n = factor_cgc_right_dict.pop("N")
            right_tmp_dict = {}
            for right_key, factor_cgc_right in factor_cgc_right_dict.items():
                right_tmp = sp.sign(factor_cgc_right) * sp.sqrt(abs(factor_cgc_right))
                right_tmp_dict[right_key] = right_tmp
            factor_cgc_right_dict["N"] = right_n
            left_right_sqrt_n = sp.sqrt(left_n * right_n)

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
                    σ_in_element = data_σ_μ.in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1, m_σ_right - 1]  # m-1得到py序
                    μ_in_element = data_σ_μ.in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1, m_μ_right - 1]
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += in_element_sum * left_tmp * right_tmp
        matrix_element = matrix_element / left_right_sqrt_n

        return True, matrix_element

    def _calc_isf_matrix(self, row_index_tmp_list, ν_st, yt_st_num_dict, data_σ_μ):
        """计算ISF的本征矩阵"""
        # TODO 根据时间反馈决定要不要上多线程/进程
        matrix_div = len(row_index_tmp_list)
        isf_matrix = sp.zeros(matrix_div)

        # 构建主对角元
        for rc, (σ_st, μ_st, β_st) in enumerate(row_index_tmp_list):
            cgc_st_tuple = (σ_st, μ_st, ν_st, β_st, 1)
            flag, isf_matrix_element = self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple, yt_st_num_dict,
                                                                     data_σ_μ)
            if not flag:
                err_msg = "get isf_matrix_element by cgc_st_tuple={}, cgc_st_tuple={}, yt_st_num_dict={}, " \
                          "data_σ_μ={}, with msg={}".format(cgc_st_tuple, cgc_st_tuple, yt_st_num_dict,
                                                            data_σ_μ, isf_matrix_element)
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
                                                                         yt_st_num_dict, data_σ_μ)
                if not flag:
                    err_msg = "get isf_matrix_element by cgc_st_tuple_left={}, cgc_st_tuple_right={}, " \
                              "yt_st_num_dict={}, data_σ_μ={} " \
                              "with msg={}".format(cgc_st_tuple_left, cgc_st_tuple_right, yt_st_num_dict,
                                                   data_σ_μ, isf_matrix_element)
                    logger.error(err_msg)
                    return False, err_msg
                isf_matrix[row, col] = isf_matrix_element
                isf_matrix[col, row] = isf_matrix_element  # 因为ISF的本征矩阵共轭

        return True, isf_matrix


class CGCHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(CGCHelper, self).__init__()
        
    def calc_cgc_dict_part_and_save_by_isf(self, isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
        # TODO 一张ISF表是可以算出对应m的完整CGC的，所以可以在里面检查N，也有可能分离save与calc
        # 3-303

        # 先按照列计算吧，后面根据情况看看是否能优化
        σ_μ_β_all_st_tuple_list = isf_square_dict["rows"]
        ν_β_list = isf_square_dict["cols"]
        isf_square_matrix = isf_square_dict["isf"]
        # isf_square_matrix_t = isf_square_matrix.T

        # TODO 待优化: 大概就是把某层最轻的循环，先包装成字典，避免重复
        for i, ν_β in enumerate(ν_β_list):
            isf_square_vector = isf_square_matrix.col(i)  # sympy库自己就是这么取行/列的
            ν, β = ν_β if isinstance(ν_β, tuple) else (ν_β, None)

            max_m_ν_st = data_st.yt_num_dict[tuple(ν_st)]
            offset_of_m_ν = self._calc_m_with_m_st(ν_st, 0, data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)
            for m_ν_st in range(1, max_m_ν_st + 1):
                single_cgc_part_start_time = time.time()

                m_ν = offset_of_m_ν + m_ν_st
                flag, cgc_square_part_dict = load_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν)
                if not flag:
                    err_msg = "get cgc_square_part_dict fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, β={}, m_ν={}, with msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ,
                                                                       ν, β, m_ν, cgc_square_part_dict)
                    logger.error(err_msg)
                    return False, err_msg
                if cgc_square_part_dict is False:  # 说明这是一个全新的
                    cgc_square_part_dict = {"N": 1}

                for σ_μ_β_all_st, isf_square in zip(σ_μ_β_all_st_tuple_list, isf_square_vector):  # TODO 这个循环里的load_cgc肯定重复了，它和ν_β无关
                    σ_st, μ_st, β_st = σ_μ_β_all_st if len(σ_μ_β_all_st) == 3 else (*σ_μ_β_all_st, None)
                    flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
                                                        is_flag_true_if_not_s_n=False)
                    if not flag:
                        err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, " \
                                  "m_ν_st={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st,
                                                             cgc_st_square_dict)
                        logger.error(err_msg)
                        return False, err_msg
                    cgc_st_square_n = cgc_st_square_dict.pop("N")

                    offset_of_m_σ = self._calc_m_with_m_st(σ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
                                                           data_st.yt_num_dict)
                    offset_of_m_μ = self._calc_m_with_m_st(μ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
                                                           data_st.yt_num_dict)
                    for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
                        cgc_key = (offset_of_m_σ + m_σ_st, offset_of_m_μ + m_μ_st)

                        if cgc_key not in cgc_square_part_dict:
                            single_cgc_square = isf_square * cgc_st_square / cgc_st_square_n
                            if single_cgc_square != 0:
                                cgc_square_part_dict[cgc_key] = isf_square * cgc_st_square / cgc_st_square_n
                        else:
                            cgc_new_part = sp.sign(isf_square) * sp.sign(cgc_st_square) \
                                           * sp.sqrt(abs(isf_square * cgc_st_square / cgc_st_square_n))
                            cgc_old_part = sp.sign(cgc_square_part_dict[cgc_key]) \
                                           * sp.sqrt(abs(cgc_square_part_dict[cgc_key]))
                            update_cgc_square = sp.sign(cgc_new_part + cgc_old_part) \
                                                * (cgc_new_part + cgc_old_part)**2
                            if update_cgc_square != 0:
                                cgc_square_part_dict[cgc_key] = update_cgc_square  # 覆盖
                            else:
                                cgc_square_part_dict.pop(cgc_key)

                # 部分的算好了，开始计算当前部分的N
                cgc_square_part_dict.pop("N")
                update_cgc_square_n = sum(abs(i) for i in cgc_square_part_dict.values())
                cgc_square_part_dict["N"] = update_cgc_square_n

                single_cgc_part_speed_time = int(time.time() - single_cgc_part_start_time)
                flag, msg = save_cgc(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_part_dict,
                                     single_cgc_part_speed_time)
                if not flag:
                    err_msg = "save_cgc fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, β={}, m_ν={}, cgc_square_part_dict={} with " \
                              "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_part_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg
                # logger.debug("ν_st={}, cgc_square_part_dict={}".format(ν_st, cgc_square_part_dict))

        return True, None
