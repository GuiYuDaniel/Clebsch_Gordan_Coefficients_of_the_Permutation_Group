# -*- coding: utf-8 -*-
"""
this code for creating the meta and the symmetry Young Diagram combination
计算Sn的元组合和对称组合
元组合指的是σμν按照杨图顺序循环，所有第一个无法被已有组合通过对称实现的组合
对称组合指的是通过对称性作用于元组合，能得到的那些组合
"""


import json
import time
from itertools import product, combinations_with_replacement, compress
from conf.cgc_config import default_s_n, group_d3, group_k4, min_s_n_of_sym
from core.young_diagrams import load_young_diagrams, calc_young_diagram_tilde
from core.cg_series import load_cg_series
from core.cgc_utils.cgc_db_typing import SYMInfo
from core.cgc_utils.cgc_local_db import get_meta_σμν_file_name, get_sym_σμν_file_name
from core.cgc_utils.cgc_local_db import get_symmetry_combination_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


def create_symmetry_combination(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储元组合及对称组合实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_symmetry_combination get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_symmetry_combination_finish_s_n()
    if not flag:
        err_msg = "get_symmetry_combination_finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} symmetry combination had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc symmetry combination s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    ϵ_key2groups_dict, groups2ϵ_key_dict = calc_ϵ_map_dicts()

    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        flag, meta_σμν_tuple_list = calc_symmetry_combination(s_i, groups2ϵ_key_dict)
        if not flag:
            err_msg = "calc_symmetry_combination fail by s_n={}, groups2ϵ_key_dict={} with " \
                      "msg={}".format(s_i, groups2ϵ_key_dict, msg)
            logger.error(err_msg)
            return False, err_msg
        flag, msg = save_meta_σμν(s_i, meta_σμν_tuple_list)
        if not flag:
            err_msg = "save_meta_σμν fail by s_n={}, meta_σμν_tuple_list={} with " \
                      "msg={}".format(s_i, meta_σμν_tuple_list, msg)
            logger.error(err_msg)
            return False, err_msg

        # 别忘了也要更新Finish_Sn
        s_i_speed_time = int(time.time() - s_i_start_time)
        flag, msg = save_symmetry_combination_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_symmetry_combination_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_symmetry_combination s_n from {} to {} done, return True, "
                "finish_s_n={}, using time={}s".format(finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def calc_symmetry_combination(s_n, groups2ϵ_key_dict):
    """
    计算所有σμν组合中的meta_σμν，同时把sym_σμν那部分归类，并记录meta_key
    """
    # 准备
    # 准备yd_list
    _, yd_list = load_young_diagrams(s_n, is_flag_true_if_not_s_n=False)
    # 准备yd的共轭yd_tilde
    yd_tilde_list = []
    for yd in yd_list:
        _, yd_tilde = calc_young_diagram_tilde(yd, is_check_yd=False)
        yd_tilde_list.append(yd_tilde)

    def get_tilde(x):
        return yd_tilde_list[yd_list.index(x)]
    # get_tilde = lambda x: yd_tilde_list[yd_list.index(x)]
    '''PEP:do not assign a lambda expression, use a def 但是，lambda不是比def开销小么，为什么不建议这种使用了呢'''

    # 开始
    meta_σμν_tuple_list = []
    all_j_σμν_set_tmp = set()  # 所有出现的组合都装这里也不是很大
    for meta_σ, meta_μ in combinations_with_replacement(yd_list, 2):  # 1/3，不是元组合的μσ不进入循环体
        _, cg_series_list = load_cg_series(s_n, meta_σ, meta_μ)
        for meta_ν in compress(yd_list, cg_series_list):  # 2/3，不是CG组合的σμν也必然不是元组合
            meta_σμν = (meta_σ, meta_μ, meta_ν)
            json_meta_σμν = tuple(json.dumps(yd) for yd in meta_σμν)
            if json_meta_σμν in all_j_σμν_set_tmp:  # 3/3，这个组合出现过，说明不是元组合
                continue
            meta_σμν_tuple_list.append(meta_σμν)  # 还是按顺序记录元组合
            # 对单个meta_σμν求sym_σμν_dict
            json_sym_σμν_list_tmp = []  # 省去每轮现做一次list(json_sym_σμν_dict.keys())
            json_sym_σμν_dict = {}
            for d3, k4 in product(group_d3, group_k4):
                sym_σμν = tuple(meta_σμν[d] if k is False else get_tilde(meta_σμν[d]) for d, k in zip(d3, k4))
                json_sym_σμν = tuple(json.dumps(yd) for yd in sym_σμν)
                meta_key = groups2ϵ_key_dict[(d3, k4)]  # 表示meta_σμν经过meta_key变换为sym_σμν
                if json_sym_σμν in json_sym_σμν_list_tmp:
                    # 出现过，只需在已有的key-value下，增补str来源
                    json_sym_σμν_dict[json_sym_σμν].append(meta_key)
                else:
                    # 未出现过，就需要新增key-value
                    json_sym_σμν_dict[json_sym_σμν] = [meta_key]
                    # 还需要增补进json_sym_σμν_list_tmp和all_j_σμν_set
                    json_sym_σμν_list_tmp.append(json_sym_σμν)
                    all_j_σμν_set_tmp |= {json_sym_σμν}  # 怕漏不怕重，所以使用并集
            # 做完单个meta_σμν的sym_σμν_dict，保存
            flag, msg = save_single_sym_σμν(s_n, meta_σ, meta_μ, meta_ν, json_sym_σμν_dict)
            if not flag:
                err_msg = "save_single_sym_σμν fail by s_n={}, σ={}, μ={}, ν={}, json_sym_σμν_dict={} with " \
                          "msg={}".format(s_n, meta_σ, meta_μ, meta_ν, json_sym_σμν_dict, msg)
                logger.error(err_msg)
                return False, err_msg

    return True, meta_σμν_tuple_list


def save_meta_σμν(s_n: int, meta_σμν_tuple_list: list):
    """
    meta_σμν的落盘格式为：
    <CG>/symmetry_info/Sn/meta_σμν.pkl  ->
    {
    "file_name": "Sn/meta_σμν",
    "data": meta_σμν_tuple_list,
    "flags": {}
    }
    meta_σμν_tuple_list = [(meta_σ, meta_μ, meta_ν), ...]
    例如：
    # <CG>/symmetry_info/S3/meta_σμν.pkl
    # {"data": [([3], [3], [3]), ([3], [2, 1], [2, 1]), ([2, 1], [2, 1], [2, 1])]}
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(meta_σμν_tuple_list, list):
        err_msg = "sym_σμν_dict={} with type={} must be list".format(meta_σμν_tuple_list, type(meta_σμν_tuple_list))
        logger.error(err_msg)
        return False, err_msg

    db_info = SYMInfo(s_n)
    _, file_name = get_meta_σμν_file_name(s_n)

    table = {"file_name": file_name,
             "data": meta_σμν_tuple_list,
             "flags": {}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_single_sym_σμν(s_n: int, σ: list, μ: list, ν: list, sym_σμν_dict: dict):
    """
    sym_σμν_dict的落盘格式为：

    <CG>/symmetry_info/Sn/[σ][μ][ν]_symmetries.pkl  ->
    {
    "file_name": "Sn/[σ][μ][ν]_symmetries",
    "data": sym_σμν_dict,
    "flags": {}
    }
    sym_σμν_dict = {(j_σ_s, j_μ_s, j_ν_s): [meta_key]}  # 来源可以不唯一
    例如：
    # <CG>/symmetry_info/S3/[2, 1][2, 1][2, 1]_symmetries.pkl
    # {"data": {('[2, 1]', '[2, 1]', '[2, 1]'): ["σμν", "μσν", ...]}}
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν, type(σ), type(μ), type(ν))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(sym_σμν_dict, dict):
        err_msg = "sym_σμν_dict={} with type={} must be dict".format(sym_σμν_dict, type(sym_σμν_dict))
        logger.error(err_msg)
        return False, err_msg

    db_info = SYMInfo(s_n)
    _, file_name = get_sym_σμν_file_name(s_n, σ, μ, ν)

    table = {"file_name": file_name,
             "data": sym_σμν_dict,
             "flags": {}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_symmetry_combination_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_symmetry_combination_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = SYMInfo(0)
    _, finish_file_name = get_symmetry_combination_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {
                           "S{}".format(s_n): s_n_speed_time
                       }}}
    if finish_s_n_before == min_s_n_of_sym - 1:
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


def get_symmetry_combination_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n=0
    """
    _, finish_file_name = get_symmetry_combination_finish_s_n_name()
    flag, data = SYMInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get_symmetry_combination_finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_sym - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_sym:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_sym, data)
        return False, err_msg


def load_meta_σμν(s_n: int, is_flag_true_if_not_s_n=True):
    """
    取得s_n下meta_σμν_tuple_list
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_meta_σμν_file_name(s_n)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = SYMInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query meta_σμν with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        meta_σμν_tuple_list = data.get("data")
        if meta_σμν_tuple_list and isinstance(meta_σμν_tuple_list, list):
            # 只检查有没有 不对内容做检查了
            return True, meta_σμν_tuple_list  # bingo！

        else:
            err_msg = "meta_σμν_tuple_list queried from db, but cannot get meta_σμν_tuple_list from data" \
                      "with data={}, meta_σμν_tuple_list={} from db".format(data, meta_σμν_tuple_list)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist meta_σμν_tuple_list db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def load_sym_σμν(s_n: int, σ: list, μ: list, ν: list, is_flag_true_if_not_s_n=True):
    """
    取得s_n下sym_σμν_dict
    如果没有，根据is_flag_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_sym:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_sym)
        logger.error(err_msg)
        return False, err_msg
    if not all(isinstance(yd, list) for yd in [σ, μ, ν]):
        err_msg = "all [σ={}, μ={}, ν_st={}] must be list but type [{}, {}, {}]".format(
            σ, μ, ν, type(σ), type(μ), type(ν))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_sym_σμν_file_name(s_n, σ, μ, ν)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = SYMInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query sym_σμν with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        sym_σμν_dict = data.get("data")
        if sym_σμν_dict and isinstance(sym_σμν_dict, dict):
            # 只检查有没有 不对内容做检查了
            return True, sym_σμν_dict  # bingo！

        else:
            err_msg = "sym_σμν_dict queried from db, but cannot get sym_σμν_dict from data" \
                      "with data={}, sym_σμν_dict={} from db".format(data, sym_σμν_dict)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist sym_σμν_dict db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def spell_ϵ_key(d3, k4, default_meta="σμν", default_sym=""):
    """根据d3, k4拼ϵ_key"""
    for d, k in zip(d3, k4):
        default_sym += "{}{}".format(default_meta[d], "~" if k is True else "")
    return default_sym


def calc_ϵ_map_dicts():
    """计算全部24种ϵ的str表示与数学表示间的映射字典"""
    ϵ_key2groups_dict = {}  # {"σ~μ~ν": ((0, 1, 2), (True, True, False))}
    groups2ϵ_key_dict = {}  # {((2, 0, 1), (False, True, True)): "νσ~μ~"}
    for d3, k4 in product(group_d3, group_k4):
        ϵ_key = spell_ϵ_key(d3, k4)
        ϵ_key2groups_dict[ϵ_key] = (d3, k4)
        groups2ϵ_key_dict[(d3, k4)] = ϵ_key
    return ϵ_key2groups_dict, groups2ϵ_key_dict
