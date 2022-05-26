# -*- coding: utf-8 -*-
"""
this code for creating CG series by Sn, young_diagrams, kai
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第10节：置换群的内积和CG序列
# 通常不同的自由度有各自的置换群。当组合起来的时候（例如q=(x, s)），就是把分量线性组合。
# CG序列就是线性组合的系数

# 公式
# [ν_1] * [ν_2] = sum_ν{(ν_1 ν_2 ν)[v]}
# (ν_1 ν_2 ν) = sum_i{g_i χ_i^ν_1 χ_i^ν_2 χ_i^ν} / g
# 其中，
# Sn(x)的IR[ν_1] 和 Sn(s)的IR[ν_2] 的乘积表示中包含 Sn(q)的IR[ν] （CG序列正是它的线性组合系数）


import copy
import math
import numpy as np
import time
from functools import lru_cache
from conf.cgc_config import default_s_n
from core.characters_and_gi import load_characters_and_gi
from core.young_diagrams import load_young_diagrams
from core.cgc_utils.cgc_db_typing import CGSeriesInfo
from core.cgc_utils.cgc_local_db import get_cg_series_file_name, get_cg_series_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


@lru_cache(maxsize=32)
def _calc_factorial_n(n):
    return math.factorial(n)


def create_cg_series(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储CG序列实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int):
        err_msg = "s_n={} with type={} must be int".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_cg_series get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_cg_series_finish_s_n()
    if not flag:
        err_msg = "get cg_series finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} cg_series had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc cg_series s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的CG series并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()

        # 数据准备
        # 这个yd_list是作为：1，matrix行和列的序号引入的；2，是作为CG序列的三个构型引入的
        flag, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db for create_cg_series meet error with s_i={}, msg={}".format(s_i, yd_list)
            logger.error(err_msg)
            return False, err_msg
        flag, characters_and_gi_dict = load_characters_and_gi(s_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get characters_and_gi_dict for create_cg_series meet error with s_i={}, msg={}".format(
                s_i, characters_and_gi_dict)
            logger.error(err_msg)
            return False, err_msg
        characters_matrix = characters_and_gi_dict.get("character")
        if not isinstance(characters_matrix, np.ndarray) or characters_matrix.size != len(yd_list) * len(yd_list):
            err_msg = "characters_matrix={} with type={} must be np.ndarray and size={} must be {} but not".format(
                characters_matrix, type(characters_matrix), characters_matrix.size, len(yd_list) * len(yd_list))
            logger.error(err_msg)
            return False, err_msg
        gi_array = characters_and_gi_dict.get("gi")
        if not isinstance(gi_array, np.ndarray) or gi_array.size != len(yd_list):
            err_msg = "gi_array={} with type={} must be np.ndarray and size={} must be {} but not".format(
                gi_array, type(gi_array), gi_array.size, len(yd_list))
            logger.error(err_msg)
            return False, err_msg
        factorial_n = _calc_factorial_n(s_i)  # factorial_n根据物理意义，一定会等于sum(gi)

        '''
        这里留的是按照物理公式，直译出的代码，没有下面优化过的快（但是全局来看，这个耗时量级，也不是很慢）
        # 1 # 30 * 30 div 0.1608600616455078
        # 2 # 30 * 30 div 0.018501758575439453
        # 结论 还是快10倍的，就是没啥用
        # 1，源码
        order_list = []
        for i in yd_list:  # [σ]
            kai_i = characters_matrix[yd_list.index(i), :]
            for j in yd_list:  # [μ]
                kai_j = characters_matrix[yd_list.index(j), :]
                order_list = []
                for k in yd_list:  # [ν]
                    kai_k = characters_matrix[yd_list.index(k), :]
                    order_k = np.dot(gi_array, (kai_i * kai_j * kai_k)) / factorial_n
                    order_k = int(order_k)
                    order_list.append(order_k)
                order_array = np.array(order_list)
                SAVE_THIS_DATA()
        '''

        # TODO 更进一步的张量优化，可以作为最佳实践去思考
        for yd_left in yd_list:  # [σ]
            kai_left = characters_matrix[yd_list.index(yd_left), :]
            for yd_right in yd_list:  # [μ]
                single_cg_series_start_time = time.time()
                kai_right = characters_matrix[yd_list.index(yd_right), :]

                k3_array = kai_left * kai_right * characters_matrix
                order_array = np.sum(gi_array * k3_array, axis=1) / factorial_n  # 按照行求和 # 按照[ν]的Yamanouchi序排列
                order_array = order_array.astype(np.int)

                # 既然没问题了，那就存（别忘了也要更新Finish_Sn）
                single_cg_series_speed_time = int(time.time() - single_cg_series_start_time)
                flag, msg = save_cg_series(s_i, yd_left, yd_right, order_array, single_cg_series_speed_time)
                if not flag:
                    err_msg = "save CG series meet error with s_i={}, yd_left={}, yd_right={}, order_array={}, " \
                              "msg={}".format(s_i, yd_left, yd_right, order_array, msg)
                    logger.error(err_msg)
                    return False, err_msg
        s_i_speed_time = int(time.time() - s_i_start_time)
        flag, msg = save_cg_series_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_cg_series_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_cg_series s_n from {} to {} done, return True, finish_s_n={}, "
                "using time={}s".format(finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_cg_series(s_n: int, yd_1: list, yd_2: list, cg_series: np.ndarray, speed_time: int):
    """
    CG序列的落盘格式为：
    <CG>/cg_series_info/Sn/[σ]_[μ].pkl        ->
    {
    "file_name": "Sn/[σ]_[μ]",
    "data": {"[ν]": cg_series},
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；
    [ν]表示内积后的可能构型；
    cg_series是[ν]的线性组合系数；
    speed_time表示计算用时（秒）

    例如：
    <CG>/cg_series_info/Sn/[3]_[2, 1].pkl
    np.array([0, 1, 0])
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd_1, list) or not isinstance(yd_2, list):
        err_msg = "yd_1={} and yd_2={} with type={}, {} all must be list".format(yd_1, yd_2, type(yd_1), type(yd_2))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(cg_series, np.ndarray):
        err_msg = "get yd_list={} with type={} must be np.ndarray".format(cg_series, type(cg_series))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = CGSeriesInfo(s_n)
    _, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
    table = {"file_name": file_name,
             "data": cg_series,
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_cg_series_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_cg_series_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = CGSeriesInfo(0)
    _, finish_file_name = get_cg_series_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": np.array([0]),
             "flags": {"finish_s_n": s_n,
                       "history_times": {"S{}".format(s_n): s_n_speed_time},
                       "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}}

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


def get_cg_series_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 1
    """
    _, finish_file_name = get_cg_series_finish_s_n_name()
    flag, data = CGSeriesInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get cg_series finish_s_n meet error with finish_file_name={}".format(finish_file_name)
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


def load_cg_series(s_n: int, yd_1: list, yd_2: list, is_flag_true_if_not_s_n=True):
    """
    取得s_n下[σ]*[μ]的CG序列（按照[ν]排列的array）
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd_1, list) or not isinstance(yd_2, list):
        err_msg = "yd_1={} and yd_2={} with type={}, {} all must be list".format(yd_1, yd_2, type(yd_1), type(yd_2))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_cg_series_file_name(s_n, yd_1, yd_2)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, yd_1={}, yd_2={} because {}".format(s_n, yd_1, yd_2, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = CGSeriesInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query cg_series with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        cg_series_array = data.get("data")
        if isinstance(cg_series_array, np.ndarray) and cg_series_array.size > 0:
            # 只检查有没有 不对内容做检查了
            return True, cg_series_array  # bingo！

        else:
            err_msg = "cg_series_array queried from db, but cannot get it from data" \
                      "with data={}, cg_series_array={} from db".format(data, cg_series_array)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist cg_series_array db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def calc_single_cg_series():
    """
    本例整体计算，没有single
    """
    pass


"""
备选db格式：
1，Sn/[σ]_[μ]/[ν]: int
2，Sn/[σ]_[μ]: dict([ν]: int) [当前！]
3，Sn: np.ndarray 3维依次是 [σ] [μ] [ν]，按照Yamanouchi序排列矩阵
"""
