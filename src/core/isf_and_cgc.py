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
import numpy as np
from functools import lru_cache, partial
from itertools import product, islice
from conf.cgc_config import default_s_n, isf_0_error_value
from core.young_diagrams import load_young_diagrams, is_young_diagram
from core.young_tableaux import load_young_table_num
from core.branching_laws import load_branching_law
from core.yamanouchi_matrix import load_yamanouchi_matrix
from core.cg_series import load_cg_series
from core.eigenvalues import load_eigenvalues
from core.cgc_utils.cgc_db_typing import ISFInfo, CGCInfo
from core.cgc_utils.cgc_local_db import get_isf_file_name, get_isf_finish_s_n_name
from core.cgc_utils.cgc_local_db import get_cgc_file_name, get_cgc_finish_s_n_name
from utils.log import get_logger
# φ


logger = get_logger(__name__)


def create_isf_and_cgc(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储isf和cgc实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    # TODO restart
    if not isinstance(s_n, int):
        err_msg = "s_n={} with type={} must be int".format(s_n, type(s_n))
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
    if isf_finish_s_n != cgc_finish_s_n:
        err_msg = "util restart code enable, isf_finish_s_n={} must eq cgc_finish_s_n={}".format(isf_finish_s_n,
                                                                                                 cgc_finish_s_n)
        logger.error(err_msg)
        return False, err_msg
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
        isf_func.enable_now_s_n(s_i)
        cgc_func.enable_now_s_n(s_i)

        # Sn循环可以得到的数据(如果可以从上一次循环继承，则继承)
        data_st_st_yt_num_dict = data_st.yt_num_dict if data_st is not None else {}
        data_st = data_si if data_si is not None else DataHelper(s_i - 1)
        data_si = DataHelper(s_i)

        # σ μ
        for σ, μ in product(data_si.yd_list, repeat=2):  # [σ], [μ]双循环
            # TODO σ, μ = μ, σ的情况，看看是真算用来消减误差，还是不算用来节约算力
            # σ_μ_start_time = time.time()

            # σ, μ循环可以得到的数据
            data_σ_μ = ΣMDataHelper(s_i, σ, μ, data_si, data_st)

            for ν_st in data_st.yd_list:  # [ν']循环
                # # ISF
                single_isf_start_time = time.time()
                # row_index_list = isf_func.calc_row_indexes(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                #                                                data_σ_μ.cg_series_st_list_dict, ν_st, data_st.yd_list)
                # isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                #                    "cols":[[ν], ([ν], β), ...],  # 有自由度tuple，无自由度list
                #                    "isf": isf_square_matrix}  # np.array([len(rows), len(cols)], dtype=float)
                flag, isf_square_dict = isf_func.calc_isf_dict(ν_st, data_si, data_st, data_σ_μ,
                                                               data_st_st_yt_num_dict)
                if not flag:
                    err_msg = "calc isf_square_dict meet error by ν_st={}, data_si={}, data_st={}, data_σ_μ={}, " \
                              "data_st_st_yt_num_dict={} with msg={}".format(ν_st, data_si, data_st, data_σ_μ, 
                                                                             data_st_st_yt_num_dict, isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
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
                flag, msg = cgc_func.calc_cgc_dict_part_and_save(isf_square_dict, ν_st, data_si, data_st, data_σ_μ)
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
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
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


def save_cgc_by_part(s_n: int, σ: list, μ: list, ν: list, β: (int, None), m: int, cgc_square_dict: dict,
                     speed_time: int):
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
        err_msg = "save_cgc_by_part meet error with file_name={} msg={}".format(file_name, old_data)
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
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
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

    if finish_s_n_before == 0:
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

    if finish_s_n_before == 0:
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
        # logger.debug("find no finish_s_n, return 0")
        return True, 0
    finish_s_n = data.get("flags", {}).get("finish_s_n")

    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= 1:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and > 0, with data={}".format(finish_s_n, data)
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
        return True, 0
    finish_s_n = data.get("flags", {}).get("finish_s_n")

    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= 1:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and > 0, with data={}".format(finish_s_n, data)
        return False, err_msg


def load_isf(s_n: int, σ: list, μ: list, ν_st: list, is_flag_true_if_not_s_n=True,
             output_mode="", ex_params=None):
    """支持两种load，分别是：

    1，只输入必要参数。返回包含isf矩阵及行、列index的字典；
    2.1，output_mode='single_row'。返回包含单行isf及列index
    2.2，output_mode='single_col'。返回包含单列isf及行index
    2.3，output_mode='single_isf'。返回单独一个isf
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
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
                if ex_params in isf_square_dict.get("rows", []):
                    row_index = isf_square_dict.get("rows").index(ex_params)
                    rst_dict = {"cols": isf_square_dict["cols"],
                                "single_row": isf_square_dict["isf"][row_index, :]}
                    return True, rst_dict  # bingo(1/4)！
                else:
                    err_msg = "ex_params={} should in rows={} with isf_square_dict={} but not, pls check".format(
                        ex_params, isf_square_dict.get("rows"), isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            elif output_mode == "single_col":
                if ex_params in isf_square_dict.get("cols", []):
                    cols_index = isf_square_dict.get("cols").index(ex_params)
                    rst_dict = {"rows": isf_square_dict["rows"],
                                "single_col": isf_square_dict["isf"][:, cols_index]}
                    return True, rst_dict  # bingo(2/4)！
                else:
                    err_msg = "ex_params={} should in cols={} with isf_square_dict={} but not, pls check".format(
                        ex_params, isf_square_dict.get("cols"), isf_square_dict)
                    logger.error(err_msg)
                    return False, err_msg
            elif output_mode == "single_isf":
                if ex_params[0] in isf_square_dict.get("rows", []) and ex_params[1] in isf_square_dict.get("cols", []):
                    row_index = isf_square_dict.get("rows").index(ex_params)
                    cols_index = isf_square_dict.get("cols").index(ex_params)
                    rst = isf_square_dict["isf"][row_index, cols_index]
                    return True, rst  # bingo(3/4)！
                else:
                    err_msg = "ex_params[0]={} should in rows={} and ex_params[1]={} should in cols={} " \
                              "with isf_square_dict={} but not, pls check".format(
                        ex_params[0], isf_square_dict.get("rows"),
                        ex_params[1], isf_square_dict.get("cols"), isf_square_dict)
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
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
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
        if not isinstance(s_n, int) or s_n < 2:
            raise Exception("s_n={} must be int and >= 2".format(s_n))
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
        pass

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
        if not isinstance(s_k, int) or s_k < 2:
            raise Exception("s_k={} must be int and >= 2".format(s_k))
        self.s_n = s_k  # 既可以实例化Sn，也可以实例化St

        self.yd_list = None
        self.bl_yd_list_dict = None
        self.yt_num_dict = None
        self.eigenvalue_list = None
        self._init_data()

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
        if not isinstance(s_n, int) or s_n < 2:
            raise Exception("s_n={} must be int and >= 2".format(s_n))
        self.s_n = s_n
        self.s_t = s_n - 1

    @staticmethod
    def _calc_m_with_m_st(yd_st, m_st, bl_yd_list, yt_st_num_dict):
        """通过m'计算m，这回，我们可以使用分支律以及yt_num了; 也可以通过m_st=0，计算偏移量"""
        # TODO 和core.young_tableaux.py: quickly_calc_young_table_in_decreasing_page_order比较一下
        bl_yd_index = bl_yd_list.index(yd_st)  # 必须成立
        # 对yd_st的前所有分支的yt_st_num以及本项m_st求和
        m = sum([yt_st_num_dict[tuple(bl_yd)] for bl_yd in bl_yd_list[:bl_yd_index]], m_st)
        return m


class ISFHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(ISFHelper, self).__init__()

    def calc_isf_dict(self, ν_st, data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict):
        """计算未经相位调整的ISF
        
        正则ISF索引的全部参数为：σ σ' μ μ' ν β ν' β'
        表示：|σ σ'> * |μ μ'> 的结果中，|νβ ν'β'>，的ISF系数平方"""
        # 只象征性检查Sn、St
        if not all(self.s_n == s_n for s_n in [data_sn.s_n, data_σ_μ.s_n]) \
                or not all(self.s_t == s_t for s_t in [data_st.s_n, data_σ_μ.s_t]):
            err_msg = "input wrong data with s_n={} should eq all(data_sn.s_n={}, data_σ_μ.s_n={}) " \
                      "and s_t={} should eq all(data_st.s_n={}, data_σ_μ.s_t={}) but not".format(
                self.s_n, data_sn.s_n, data_σ_μ.s_n, self.s_t, data_st.s_n, data_σ_μ.s_t)
            logger.error(err_msg)
            return False, err_msg

        row_index_tmp_list = self.calc_row_indexes_tmp(data_σ_μ.bl_yds_of_σ, data_σ_μ.bl_yds_of_μ,
                                                       data_σ_μ.cg_series_st_list_dict, ν_st, data_st.yd_list)

        isf_matrix = self._calc_isf_matrix(row_index_tmp_list, ν_st,
                                           data_st.bl_yd_list_dict, data_st.yt_num_dict,
                                           data_σ_μ.in_matrix_σ_dict, data_σ_μ.in_matrix_μ_dict)
        eigenvalues, eigenvectors = np.linalg.eigh(isf_matrix)  # eigh适用于复共轭情况
        eigenvalues_int = [int(np.around(i, 0)) for i in eigenvalues]  # 理论上这里的本征值必为整数  # TODO 要检查
        # 把本征值和本征矢量改为降序，T是为了适应for循环逻辑。for循环是按照行来取，而np.linalg.eigh按列取
        # 既，np.linalg.eigh算出来的eigenvectors需要eigenvectors[:, index(eigenvalues)]取值，与for循环逻辑不同，所以需要.T调整
        eigenvalues_int = eigenvalues_int[::-1]
        eigenvectors_t = eigenvectors[::-1].T  # 将默认先列的numpy结果转化为先行

        # 先进行施密特正交归一化，并计算出ν_list, β_list
        soe_vectors, β_list = self._calc_schmidt_orthogonalization_eigenvectors_and_β_list(eigenvalues_int,
                                                                                           eigenvectors_t)

        λ_ν_st = data_st.eigenvalue_list[data_st.yd_list.index(ν_st)]
        # 开始计算ISF
        isf_square_tmp_dict = {"ν_tmp": [],
                               "β_tmp": [],
                               "isf_tmp": []}
        for e_value, s_vector, β in zip(eigenvalues_int, soe_vectors, β_list):  # 通过e_value建立ν
            # TODO [3,1,1]_[3,1,1] [3,1]和[2,2]的细节打印出来，观察β
            # 计算本征值对应的ν
            λ_ν = e_value + λ_ν_st
            ν = self._calc_ν_by_λ_and_bl(λ_ν, data_sn.eigenvalue_list, data_sn.yd_list,
                                         data_sn.bl_yd_list_dict, ν_st)

            # 调整eigenvector相位  # TODO 看看可否通过调整e_value循环顺序进行优化(np.linalg.eigh默认给出的本征值是升序排列)
            # TODO 这里有可能需要传入isf_square_tmp_dict，看debug吧
            flag, isf_phase = self._calc_isf_phase(s_vector, ν, β, ν_st, row_index_tmp_list,
                                                   data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict)
            if not flag:
                err_msg = "calc isf_phase meet error with s_vector={}, ν={}, β={}, ν_st={}, row_index_tmp_list={}, " \
                          "data_sn={}, data_st={}, data_σ_μ={}, data_st_st_yt_num_dict={}, " \
                          "msg={}".format(s_vector, ν, β, ν_st, row_index_tmp_list, data_sn, data_st, data_σ_μ,
                                          data_st_st_yt_num_dict, isf_phase)
                logger.error(err_msg)
                return False, err_msg

            # 对确定相位来说，最重要的不是数值，而是符号  # TODO 如果0化，小心sign(0) = 0
            phase_vector = s_vector * np.sign(isf_phase)

            # 计算单列的ISF
            isf_square = np.sign(phase_vector) * (phase_vector * phase_vector)
            isf_square_tmp_dict["ν_tmp"].append(ν)
            isf_square_tmp_dict["β_tmp"].append(β)  # β按照代码逻辑，同ν必然升序，且临近
            isf_square_tmp_dict["isf_tmp"].append(isf_square)

        # 按照顺序整理矩阵
        row_index_list = [(i[0], i[1]) if i[2] is None else i for i in row_index_tmp_list]
        isf_square_dict = {"rows": row_index_list,
                           "cols": [],
                           "isf": np.zeros([len(row_index_list), len(row_index_list)])}
        if len(isf_square_tmp_dict["ν_tmp"]) != len(row_index_list):  # 先检查理论上要求的方阵
            err_msg = "ν_tmp_list={} with len={} and row_index_list={} with len={} must same, pls check".format(
                isf_square_tmp_dict["ν_tmp"], len(isf_square_tmp_dict["ν_tmp"]),
                row_index_list, len(row_index_list))
            logger.error(err_msg)
            return False, err_msg

        for ν in data_sn.yd_list:
            while ν in isf_square_tmp_dict["ν_tmp"]:
                # 找数据
                tmp_index = isf_square_tmp_dict["ν_tmp"].index(ν)
                β = isf_square_tmp_dict["β_tmp"][tmp_index]
                single_col_index = ν if β is None else (ν, β,)
                isf_square = isf_square_tmp_dict["isf_tmp"][tmp_index]
                # 赋值
                isf_square_dict["cols"].append(single_col_index)
                isf_square_dict["isf"][:, tmp_index] = isf_square
                # 删tmp数据
                isf_square_tmp_dict["ν_tmp"].pop(tmp_index)
                isf_square_tmp_dict["β_tmp"].pop(tmp_index)
                isf_square_tmp_dict["isf_tmp"].pop(tmp_index)

        return True, isf_square_dict

    def _calc_isf_phase(self, s_vector, ν, β, ν_st, row_index_tmp_list,
                        data_sn, data_st, data_σ_μ, data_st_st_yt_num_dict):
        """
        1，将分支律的第一分支首个非零系数调整为正（绝对相位）
        2，非第一分支的，参考第一分支做相对调整（相对相位）
        """
        if data_sn.bl_yd_list_dict[tuple(ν)][0] == ν_st:  # ν_st击中ν的第一分支
            return True, self._get_absolute_phase(s_vector)
        else:  # 未击中情况，则相对相位需要参考它的分支律第一分支（ν_st_reference）的绝对相位
            ν_st_reference = data_sn.bl_yd_list_dict[tuple(ν)][0]
            # 因为首先按照ν_st开启循环，所以必须有结果
            flag, isf_col_reference_dict = load_isf(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st_reference,
                                                    output_mode="single_col", ex_params=(ν, β,),
                                                    is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "load_isf meet error with self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, ν_st_reference={}," \
                          "output_mode='single_col', ex_params={}, msg={}".format(
                    self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν_st_reference, (ν, β,), isf_col_reference_dict)
                logger.error(err_msg)
                return False, err_msg
            row_index_list_reference = isf_col_reference_dict["rows"]
            isf_square_list_reference = isf_col_reference_dict["single_col"]
            first_no_0_isf_reference = self._get_first_no_0_number_from_vector(isf_square_list_reference)
            first_no_0_isf_reference_row_index = \
                row_index_list_reference[isf_square_list_reference.index(first_no_0_isf_reference)]
            # 分别为：
            # 1,令非ν的第一分支ν_st的m'取1时，ν的m；
            # 2,令ν_st_st的m''取1时，ν_st_reference的m'；(ν_st_st是ν_st的第一分支，也是ν_st_reference的分支，但不一定是第一分支)
            # 3,令2中的ν_st_reference定yt后，m_ν_by_m_ν_st_reference的yt序
            m_ν_st = 1
            m_ν = self._calc_m_with_m_st(ν_st, m_ν_st, data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)
            ν_st_st = data_st.bl_yd_list_dict[tuple(ν_st)][0]  # TODO ν_st的第一分支能击中ν_st_reference的分支律的原因是
            m_ν_st_reference = self._calc_m_with_m_st(ν_st_st, m_ν_st,
                                                      data_st.bl_yd_list_dict[tuple(ν_st_reference)],
                                                      data_st_st_yt_num_dict)
            m_ν_by_m_ν_st_reference = self._calc_m_with_m_st(ν_st_reference, m_ν_st_reference,
                                                             data_sn.bl_yd_list_dict[tuple(ν)], data_st.yt_num_dict)

            flag, before_isf_reference = self._calc_before_isf_reference(
                ν, m_ν, m_ν_by_m_ν_st_reference, first_no_0_isf_reference_row_index, ν_st_reference,
                row_index_tmp_list, s_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ)
            if not flag:
                err_msg = "calc before_isf_reference meet error with ν={}, m_ν={}, m_ν_by_m_ν_st_reference={}, " \
                          "first_no_0_isf_reference_row_index={}, ν_st_reference={}, row_index_tmp_list={}, " \
                          "s_vector={}, ν_st={}, m_ν_st={}, data_sn={}, data_st={}, data_σ_μ={}, " \
                          "msg={}".format(ν, m_ν, m_ν_by_m_ν_st_reference, first_no_0_isf_reference_row_index,
                                          ν_st_reference, row_index_tmp_list, s_vector, ν_st, m_ν_st,
                                          data_sn, data_st, data_σ_μ, before_isf_reference)
                logger.error(err_msg)
                return False, err_msg

            return True, np.sign(before_isf_reference) * np.sign(first_no_0_isf_reference)

    def _calc_before_isf_reference(self, ν, m_ν, m_ν_by_m_ν_st_reference,
                                   first_no_0_isf_reference_row_index, ν_st_reference,
                                   row_index_tmp_list, s_vector, ν_st, m_ν_st, data_sn, data_st, data_σ_μ):
        """isf除了有书中介绍的两种根据St_ISF计算的递推方法外，还可以根据已有Sn的isf，推断一部分isf

        见式子 4-195c"""
        σ_st_reference, μ_st_reference = first_no_0_isf_reference_row_index[0], first_no_0_isf_reference_row_index[1]
        β_st_reference = first_no_0_isf_reference_row_index[3] if len(first_no_0_isf_reference_row_index) == 3 else None
        flag, cgc_square_st_reference_dict = load_cgc(self.s_t, σ_st_reference, μ_st_reference, ν_st_reference,
                                                      β_st_reference, m_ν_by_m_ν_st_reference,
                                                      is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get cgc_square_dict with self.s_t={}, σ_st_reference={}, μ_st_reference={}, " \
                      "ν_st_reference={}, β_st_reference={}, cgc_square_st_reference_dict={} meet error with " \
                      "msg={}".format(self.s_t, σ_st_reference, μ_st_reference, ν_st_reference, β_st_reference,
                                      m_ν_by_m_ν_st_reference, cgc_square_st_reference_dict)
            logger.error(err_msg)
            return False, err_msg
        cgc_square_reference_n = cgc_square_st_reference_dict.pop("N")

        sum_3_loop = 0
        for σ_st, μ_st, β_st, s_vector_element in zip(row_index_tmp_list, s_vector):
            flag, cgc_square_st_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, is_flag_true_if_not_s_n=False)
            if not flag:
                err_msg = "get cgc_square_st_dict with self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, m_ν_st={} " \
                          "meet error with msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν_st, cgc_square_st_dict)
                logger.error(err_msg)
                return False, err_msg
            cgc_square_dict_n = cgc_square_st_dict.pop("N")
            sum_3_loop_part = 0
            for (m_σ_st_reference, m_μ_st_reference), cgc_square_st_reference_element \
                    in cgc_square_st_reference_dict.items():
                m_σ_reference = self._calc_m_with_m_st(σ_st_reference, m_σ_st_reference,
                                                       data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)], data_st.yt_num_dict)
                m_μ_reference = self._calc_m_with_m_st(μ_st_reference, m_μ_st_reference,
                                                       data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)], data_st.yt_num_dict)
                for (m_σ_st, m_μ_st), cgc_square_st_element in cgc_square_st_dict.items():
                    m_σ = self._calc_m_with_m_st(σ_st, m_σ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
                                                 data_st.yt_num_dict)
                    m_μ = self._calc_m_with_m_st(μ_st, m_μ_st, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
                                                 data_st.yt_num_dict)
                    in_matrix_σ_element = data_σ_μ.in_matrix_σ_dict[(self.s_t, self.s_n)][m_σ_reference - 1, m_σ - 1]
                    in_matrix_μ_element = data_σ_μ.in_matrix_μ_dict[(self.s_t, self.s_n)][m_μ_reference - 1, m_μ - 1]
                    sum_3_loop_part += in_matrix_σ_element * in_matrix_μ_element \
                                       * np.sign(cgc_square_st_element) * np.sign(cgc_square_st_reference_element) \
                                       * np.sqrt(abs(cgc_square_st_element * cgc_square_st_reference_element
                                                     / (cgc_square_dict_n * cgc_square_reference_n)))
            sum_3_loop += sum_3_loop_part * s_vector_element

        flag, in_matrix_ν = load_yamanouchi_matrix(self.s_n, ν, (self.s_t, self.s_n,), mode="in")  # (Sn-1, Sn)的对换
        if not flag:
            err_msg = "get in_matrix_ν with s_n={}, ν={}, in_key={} meet error with " \
                      "msg={}".format(self.s_n, ν, (self.s_t, self.s_n,), in_matrix_ν)
            logger.error(err_msg)
            return False, err_msg
        in_matrix_ν_m_m_element = in_matrix_ν[m_ν - 1, m_ν_by_m_ν_st_reference - 1]
        reference_isf = sum_3_loop / in_matrix_ν_m_m_element
        return True, reference_isf

    def _calc_reference_isf_tmp(self):
        """TODO 看看有关reference这一段能不能独立出来成为一个函数"""
        pass

    @staticmethod
    def _get_first_no_0_number_from_vector(vector, error_value=isf_0_error_value):
        """提取vector中首个误差范围内非0的数字"""
        for number in vector:
            if abs(number) < error_value:
                pass
            return number
        return None

    def _get_absolute_phase(self, schmidt_orthogonalization_vector, error_value=isf_0_error_value):
        """判断绝对相位: 本函数的作用是提取首个不为零的vector的符号"""
        no_0_element = self._get_first_no_0_number_from_vector(schmidt_orthogonalization_vector)
        if no_0_element is None:
            msg = "find schmidt_orthogonalization_vector={} all approximately equal to 0, with isf_0_error_value={}," \
                  "pls check".format(schmidt_orthogonalization_vector, error_value)
            logger.warning(msg)
            raise Exception(msg)
        return np.sign(no_0_element)

    @staticmethod
    def _calc_schmidt_orthogonalization_eigenvectors_and_β_list(eigenvalues_int, eigenvectors_t):
        """将eigenvectors_t中eigenvalue重复的使用施密特正交归一化手续，并计算出相应的β"""
        # TODO 这里和old算法不一样，旧的算法名字叫schmidt第一分量却取的所有分量之和做的归一
        β_list = []
        e_list = []  # eigenvectors_schmidt_orthogonalization  # 列表就够了，不需要矩阵
        tmp_counter_dict = {}  # {e_value: {β_index: e_index}}  # 实时计数器，记录已经施密特正交化的那部分β
        for (e_index, e_value), e_vector in zip(enumerate(eigenvalues_int), eigenvectors_t):
            # e_value是ν的一一映射，所以count(e_value)也可以表示自由度β
            if eigenvalues_int.count(e_value) == 1:  # 无自由度
                β = None
                β_list.append(β)
                e_list.append(e_vector)
            else:  # 有自由度
                β = len(tmp_counter_dict.get(e_value, {})) + 1
                β_list.append(β)
                '''施密特正交化：
                设欧氏空间中向量a1，a2，a3线性无关，令
                b1 = a1
                b2 = a2 - <a2,b1>/||b1|| * b1
                b3 = a3 - <a3,b1>/||b1|| * b1 - <a3,b2>/||b2|| * b2'''
                if β == 1:
                    tmp_counter_dict[e_value] = {1: e_index}
                    e_list.append(e_vector)
                else:
                    # TODO 如果有必要，下面这段可以优化
                    minuend_vector = 0
                    for β_index in range(1, β):  # 这里的β_index指的是按照schmidt公式被去掉的前β-1个矢量的标号
                        cos_ab = sum(e_vector * e_list[tmp_counter_dict[e_value][β_index]])
                        # norm_a norm_b 都应该约等于1
                        # norm_ab = np.linalg.norm(e_vector) * np.linalg.norm(e_list[tmp_counter_dict[e_value][β_index]])
                        # minuend_vector += (cos_ab / norm_ab ) * tmp_counter_dict[e_value][β_index]
                        minuend_vector += cos_ab * tmp_counter_dict[e_value][β_index]
                    schmidt_vector = e_vector - minuend_vector
                    # TODO 这个归一化手续对误差的影响巨大哦！
                    schmidt_orthogonalization_vector = schmidt_vector / np.linalg.norm(schmidt_vector)
                    tmp_counter_dict[e_value][β] = e_index
                    e_list.append(schmidt_orthogonalization_vector)
        return e_list, β_list

    @staticmethod
    def _calc_schmidt_orthogonalization_with_counter(β, eigenvectors_t):
        """使用施密特正交归一化手续将拥有自由度的矢量归一化"""

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
    def calc_row_indexes_tmp(bl_yds_of_σ, bl_yds_of_μ, cg_series_st_list_dict, ν_st, yd_list):
        """计算ISF表格的行的意义，它是的bl_σ, bl_μ, β'的列表
        形如[([3,1],[3,1],None), ([3,1],[2,1,1],1), ([3,1],[2,1,1],2), ...]"""
        row_index_tmp_list = []  # [(bl_σ, bl_μ, β'), ...]
        for bl_σ, bl_μ in product(bl_yds_of_σ, bl_yds_of_μ):
            cg_series_st_list = cg_series_st_list_dict[(tuple(bl_σ), tuple(bl_μ))]  # 必须有
            single_ν_cg_series = cg_series_st_list[yd_list.index(ν_st)]
            if single_ν_cg_series == 0:
                continue
            part_rst_list = [(bl_σ, bl_μ, β_st) if single_ν_cg_series >= 2 else (bl_σ, bl_μ, None)
                             for β_st in range(1, single_ν_cg_series + 1)]  # [1, 2, ..., single_ν_cg_series]
            row_index_tmp_list += part_rst_list
        return row_index_tmp_list

    @lru_cache()
    def _load_cgc_with_m1_by_input_json(self, s_n, j_σ, j_μ, j_ν, β):
        """lru版本的load_cgc，通过将入参更改为json.dumps实现序列化/json.loads还原"""
        σ, μ, ν = map(json.loads, [j_σ, j_μ, j_ν])
        return load_cgc_with_m1(s_n, σ, μ, ν, β)

    def _load_cgc_with_m1_lru(self, s_n, σ, μ, ν, β):
        j_σ, j_μ, j_ν = map(json.dumps, [σ, μ, ν])
        return self._load_cgc_with_m1_by_input_json(s_n, j_σ, j_μ, j_ν, β)

    def _cgc_st_2_cgc_dict(self, cgc_st_tuple, bl_st_dict, yt_st_num_dict):
        """取得cg_st_dict，并且，根据分支律，将m_st转化为对应的m"""
        σ_st, μ_st, ν_st, β_st, _ = cgc_st_tuple
        rst_dict = {}
        flag, cgc_square_dict = self._load_cgc_with_m1_lru(self.s_t, σ_st, μ_st, ν_st, β_st)
        if not flag:
            err_msg = "load_cgc fail by s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={} with msg={}".format(
                self.s_t, σ_st, μ_st, ν_st, β_st, cgc_square_dict)
            logger.error(err_msg)
            return False, err_msg
        rst_n = cgc_square_dict.pop("N")
        bl_σ_st_list = bl_st_dict[tuple(σ_st)]
        bl_μ_st_list = bl_st_dict[tuple(μ_st)]
        for (m1_st, m2_st), value in cgc_square_dict.items():
            m1 = self._calc_m_with_m_st(σ_st, m1_st, bl_σ_st_list, yt_st_num_dict)
            m2 = self._calc_m_with_m_st(μ_st, m1_st, bl_μ_st_list, yt_st_num_dict)
            rst_dict[(m1, m2,)] = cgc_square_dict[(m1_st, m2_st,)]
        rst_dict["N"] = rst_n

        return True, rst_dict
    
    def _calc_isf_matrix_element(self, cgc_st_tuple_left, cgc_st_tuple_right,
                                 bl_st_dict, yt_st_num_dict, in_matrix_σ_dict, in_matrix_μ_dict):
        """计算ISF本征矩阵的矩阵元，其中主对角元可适当优化"""
        matrix_element = 0
        # 将St的cgc根据分支律上升到Sn，并整形成py序
        # 左
        factor_cgc_left_dict = self._cgc_st_2_cgc_dict(cgc_st_tuple_left, bl_st_dict, yt_st_num_dict)
        left_n = factor_cgc_left_dict.pop("N")
        left_tmp_dict = {}
        for left_key, factor_cgc_left in factor_cgc_left_dict.items():
            left_tmp = np.sign(factor_cgc_left) * np.sqrt(abs(factor_cgc_left))
            left_tmp_dict[left_key] = left_tmp
        factor_cgc_left_dict["N"] = left_n  # 还原字典，避免deepcopy
        # 右
        if cgc_st_tuple_right == cgc_st_tuple_left:  # 如果相等，则直接使用left结果，无需重复计算
            right_tmp_dict = left_tmp_dict
            left_right_sqrt_n = left_n
        else:
            factor_cgc_right_dict = self._cgc_st_2_cgc_dict(cgc_st_tuple_right, bl_st_dict, yt_st_num_dict)
            right_n = factor_cgc_right_dict.pop("N")
            right_tmp_dict = {}
            for right_key, factor_cgc_right in factor_cgc_right_dict.items():
                right_tmp = np.sign(factor_cgc_right) * np.sqrt(abs(factor_cgc_right))
                right_tmp_dict[right_key] = right_tmp
            factor_cgc_right_dict["N"] = right_n
            left_right_sqrt_n = np.sqrt(left_n * right_n)

        # 计算matrix_element
        for (m_σ_left, m_μ_left), left_tmp in left_tmp_dict.items():
            for (m_σ_right, m_μ_right), right_tmp in right_tmp_dict.items():
                in_element_sum = 0
                for i in range(1, self.s_n):  # 这个i是交换矩阵（in）中的i
                    σ_in_element = in_matrix_σ_dict[(i, self.s_n)][m_σ_left - 1][m_σ_right - 1]  # m-1得到py序
                    μ_in_element = in_matrix_μ_dict[(i, self.s_n)][m_μ_left - 1][m_μ_right - 1]
                    in_element_sum += σ_in_element * μ_in_element
                matrix_element += in_element_sum * left_tmp * right_tmp
        matrix_element = matrix_element / left_right_sqrt_n

        return matrix_element

    def _calc_isf_matrix(self, row_index_tmp_list, ν_st, bl_st_dict, yt_st_num_dict,
                         in_matrix_σ_dict, in_matrix_μ_dict):
        """计算ISF的本征矩阵"""
        # TODO 根据时间反馈决定要不要上多线程/进程
        matrix_div = len(row_index_tmp_list)
        isf_matrix = np.zeros([matrix_div, matrix_div])
        # 构建主对角元
        for rc, (σ_st, μ_st, β_st) in enumerate(row_index_tmp_list):
            cgc_st_tuple = (σ_st, μ_st, ν_st, β_st, 1)
            isf_matrix[rc][rc] = self._calc_isf_matrix_element(cgc_st_tuple, cgc_st_tuple,
                                                               bl_st_dict, yt_st_num_dict,
                                                               in_matrix_σ_dict, in_matrix_μ_dict)
        # 构建其他矩阵元素
        if matrix_div >= 2:
            for row, (σ_st_left, μ_st_left, β_st_left) in enumerate(row_index_tmp_list):
                cgc_st_tuple_left = (σ_st_left, μ_st_left, ν_st, β_st_left, 1)
                for col, (σ_st_right, μ_st_right, β_st_right) in enumerate(islice(row_index_tmp_list, row, None)):
                    cgc_st_tuple_right = (σ_st_right, μ_st_right, ν_st, β_st_right, 1)
                    isf_matrix[row][col] = self._calc_isf_matrix_element(cgc_st_tuple_left, cgc_st_tuple_right,
                                                                         bl_st_dict, yt_st_num_dict,
                                                                         in_matrix_σ_dict, in_matrix_μ_dict)
                    isf_matrix[col][row] = isf_matrix[row][col]  # 因为ISF的本征矩阵共轭

        return True, isf_matrix


class CGCHelper(CalcHelper):
    """这里定义了一些供模块内部使用的函数，并省略入参检查"""

    def __init__(self):
        super(CGCHelper, self).__init__()
        
    def calc_cgc_dict_part_and_save(self, isf_square_dict, ν_st, data_sn, data_st, data_σ_μ):
        # 先按照列计算吧，后面根据情况看看是否能优化
        σ_μ_β_all_st_tuple_list = isf_square_dict["rows"]
        ν_β_list = isf_square_dict["cols"]
        isf_square_matrix = isf_square_dict["isf"]  # 现在for它是按照行循环哦
        isf_square_matrix_t = isf_square_matrix.T

        # TODO 待优化: 大概就是把某层最轻的循环，先包装成字典，避免重复
        for ν_β, isf_square_vector in zip(ν_β_list, isf_square_matrix_t):
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
                    flag, cgc_st_square_dict = load_cgc(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν,
                                                        is_flag_true_if_not_s_n=False)
                    if not flag:
                        err_msg = "get cgc_st_square_dict fail by self.s_t={}, σ_st={}, μ_st={}, ν_st={}, β_st={}, " \
                                  "m_ν={}, msg={}".format(self.s_t, σ_st, μ_st, ν_st, β_st, m_ν, cgc_st_square_dict)
                        logger.error(err_msg)
                        return False, err_msg
                    cgc_st_square_n = cgc_st_square_dict.pop("N")

                    offset_of_m_σ = self._calc_m_with_m_st(σ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.σ)],
                                                           data_st.yt_num_dict)
                    offset_of_m_μ = self._calc_m_with_m_st(μ_st, 0, data_sn.bl_yd_list_dict[tuple(data_σ_μ.μ)],
                                                           data_st.yt_num_dict)
                    for (m_σ_st, m_μ_st), cgc_st_square in cgc_st_square_dict.items():
                        m_σ = offset_of_m_σ + m_σ_st
                        m_μ = offset_of_m_μ + m_μ_st

                        if (m_σ, m_μ) not in cgc_square_part_dict:
                            cgc_square_part_dict[(m_σ, m_μ)] = isf_square * cgc_st_square / cgc_st_square_n
                        else:
                            cgc_new_part = np.sign(isf_square) * np.sign(cgc_st_square) \
                                           * np.sqrt(abs(isf_square * cgc_st_square / cgc_st_square_n))
                            cgc_old_part = np.sign(cgc_square_part_dict[(m_σ, m_μ)]) \
                                           * np.sqrt(abs(cgc_square_part_dict[(m_σ, m_μ)]))
                            update_cgc_square = np.sign(cgc_new_part + cgc_old_part) \
                                                * (cgc_new_part + cgc_old_part)**2
                            cgc_square_part_dict[(m_σ, m_μ)] = update_cgc_square  # 覆盖

                # 部分的算好了，开始计算当前部分的N
                cgc_square_part_dict.pop("N")
                update_cgc_square_n = sum(abs(i) for i in cgc_square_part_dict.values())
                cgc_square_part_dict["N"] = update_cgc_square_n

                single_cgc_part_speed_time = int(time.time() - single_cgc_part_start_time)
                flag, msg = save_cgc_by_part(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_part_dict,
                                             single_cgc_part_speed_time)
                if not flag:
                    err_msg = "save_cgc_by_part fail by self.s_n={}, data_σ_μ.σ={}, data_σ_μ.μ={}, " \
                              "ν={}, β={}, m_ν={}, cgc_square_part_dict={} with " \
                              "msg={}".format(self.s_n, data_σ_μ.σ, data_σ_μ.μ, ν, β, m_ν, cgc_square_part_dict, msg)
                    logger.error(err_msg)
                    return False, err_msg

        return True, None
