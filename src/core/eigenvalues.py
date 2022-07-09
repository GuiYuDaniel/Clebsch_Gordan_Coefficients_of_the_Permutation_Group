# -*- coding: utf-8 -*-
"""
this code for creating eigenvalues of 2-cycle class by Sn, young_diagrams
and provide a function to calc eigenvalue2yds
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第5节：置换群的CSCO-Ⅱ
# 这里计算二循环类算符是为了标记寻找唯一的本征值集"航道"

# 公式：
# λ = 1/2 * sum{l}(ν_l * (ν_l - 2l + 1))
#   = 1/2 * sum{i}(ν_i * (ν_i - 2l - 1))
# 其中，l是杨盘图的行号（从1开始），ν_l是指第l行的格子数
# 因为python中，index从0开始，所以l = i + 1


import copy
import time
from conf.cgc_config import default_s_n, min_s_n_of_eigenvalue
from core.young_diagrams import load_young_diagrams, is_young_diagram
from core.cgc_utils.cgc_db_typing import EigenvaluesInfo
from core.cgc_utils.cgc_local_db import get_eigenvalues_file_name, get_eigenvalues_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)


def create_eigenvalues(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储eigenvalues实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_eigenvalue:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_eigenvalue)
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_eigenvalues get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_eigenvalues_finish_s_n()
    if not flag:
        err_msg = "get_eigenvalues_finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} eigenvalues had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc eigenvalues s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的eigenvalues并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()

        # 数据准备
        flag, yd_list = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db for create_eigenvalues meet error with s_i={}, msg={}".format(s_i, yd_list)
            logger.error(err_msg)
            return False, err_msg
        flag, eigenvalues_list = calc_eigenvalues_of_2_cycle_class_of_s_n(yd_list)
        if not flag:
            err_msg = "calc_eigenvalues_of_2_cycle_class_of_s_n meet error with s_i={}, yd_list={}, msg={}".format(
                s_i, yd_list, eigenvalues_list)
            logger.error(err_msg)
            return False, err_msg

        # 既然没问题了，那就存（别忘了也要更新Finish_Sn）
        s_i_speed_time = int(time.time() - s_i_start_time)
        flag, msg = save_eigenvalues(s_i, eigenvalues_list, s_i_speed_time)
        if not flag:
            err_msg = "save_eigenvalues meet error with s_i={}, eigenvalues_list={}, msg={}".format(
                s_i, eigenvalues_list, msg)
            logger.error(err_msg)
            return False, err_msg

        flag, msg = save_eigenvalues_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_eigenvalues_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_eigenvalues s_n from {} to {} done, return True, finish_s_n={}, "
                "using time={}s".format(finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_eigenvalues(s_n: int, eigenvalues_list: list, speed_time: int):
    """
    二循环类的本征值的落盘格式为：
    <CG>/eigenvalues_info/Sn.pkl         ->
    {
    "file_name": "Sn",
    "data": eigenvalues_list,             # list(int)
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    eigenvalues_list是yd按照Yamanouchi序排列对应的二循环类的本征值；
    speed_time表示计算用时（秒）

    例如：
    <CG>/eigenvalues_info/S6.pkl
    [15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15]
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_eigenvalue:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_eigenvalue)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(eigenvalues_list, list):
        err_msg = "get eigenvalues_list={} with type={} must be list".format(eigenvalues_list, type(eigenvalues_list))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = EigenvaluesInfo(s_n)
    _, file_name = get_eigenvalues_file_name(s_n)
    table = {"file_name": file_name,
             "data": eigenvalues_list,
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_eigenvalues_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_eigenvalue:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_eigenvalue)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_eigenvalues_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = EigenvaluesInfo(0)
    _, finish_file_name = get_eigenvalues_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": [],
             "flags": {"finish_s_n": s_n,
                       "history_times": {"S{}".format(s_n): s_n_speed_time},
                       "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}}

    if finish_s_n_before == min_s_n_of_eigenvalue - 1:
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


def get_eigenvalues_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n = 0
    """
    _, finish_file_name = get_eigenvalues_finish_s_n_name()
    flag, data = EigenvaluesInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get_eigenvalues_finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_eigenvalue - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_eigenvalue:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_eigenvalue, data)
        return False, err_msg


def load_eigenvalues(s_n: int, is_flag_true_if_not_s_n=True):
    """
    取得s_n下所有young diagrams的本征值列表
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_eigenvalue:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_eigenvalue)
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_eigenvalues_file_name(s_n)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = EigenvaluesInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query eigenvalues with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        eigenvalues_list = data.get("data")
        if isinstance(eigenvalues_list, list) and eigenvalues_list:
            # 只检查有没有 不对内容做检查了
            return True, eigenvalues_list  # bingo！

        else:
            err_msg = "eigenvalues_list queried from db, but cannot get it from data" \
                      "with data={}, eigenvalues_list={} from db".format(data, eigenvalues_list)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist eigenvalues_list db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def get_yds_list_by_eigenvalue(s_n: int, eigenvalue: int):
    """
    根据提供的本征值以及s_n，给出所有对应的yds_list（不总是唯一）
    """
    # TODO
    pass


def get_all_eigenvalue2yd_dict():
    """
    根据使用情况决定哪个函数更好
    """
    # TODO
    pass


def _calc_single_eigenvalue_of_2_cycle_class_of_yd(yd: list):
    """
    公式：
    λ = 1/2 * sum{l}(ν_l * (ν_l - 2l + 1))
      = 1/2 * sum{i}(ν_i * (ν_i - 2i - 1))
    其中，l是杨盘图的行号（从1开始），ν_l是指第l行的格子数
    因为python中，index从0开始，所以l = i + 1
    """
    if not is_young_diagram(yd):
        err_msg = "cannot calc eigenvalue with wrong young_diagram={}".format(yd)
        logger.error(err_msg)
        return False, err_msg
    if yd == [1]:
        return True, 1  # 注意，二循环类算符，本应从S2开始，其S1的值不是计算的，而是定义的延拓

    rst = 0
    for row in range(len(yd)):
        # rst_part = yd[row] * (yd[row] + 1 - 2 * (row + 1)) / 2
        rst_part = yd[row] * (yd[row] - 2 * row - 1)  # /2放外面
        rst += rst_part
    rst = int(rst / 2)
    return True, rst


def calc_eigenvalues_of_2_cycle_class_of_s_n(yd_list: list):
    if not isinstance(yd_list, list):
        err_msg = "yd_list={} for calc_eigenvalues_of_2_cycle_class must be list but not".format(yd_list)
        logger.error(err_msg)
        return False, err_msg

    eigenvalues_list = []
    for single_yd in yd_list:
        flag, single_eigenvalue = _calc_single_eigenvalue_of_2_cycle_class_of_yd(single_yd)
        if not flag:
            err_msg = "_calc_single_eigenvalue_of_2_cycle_class_of_yd meet error with single_yd={}, msg={}".format(
                single_yd, single_eigenvalue)
            logger.error(err_msg)
            return False, err_msg
        eigenvalues_list.append(single_eigenvalue)

    if len(yd_list) != len(eigenvalues_list):
        err_msg = "len(yd_list)={} should eq len(eigenvalues_list)={} but not with yd_list={}, eigenvalues_list={}" \
                  "".format(yd_list, eigenvalues_list, type(yd_list), type(eigenvalues_list))
        logger.error(err_msg)
        return False, err_msg

    return True, eigenvalues_list
