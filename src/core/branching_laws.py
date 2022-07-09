# -*- coding: utf-8 -*-
"""
this code for creating Branching Laws by Sn and young_diagram
计算Sn杨图的分支律
"""

# \subsection{分支律}
#
# 置换群$S_n$的$IR[\nu]$相对于其子群$S_{n-1} \equiv S_{n-1}(1,2,...,n-1)$来说，一般是可约的。
# 传统置换群表示论证明了，若用杨图$[\nu]$标志$S_n$群的不可约表示：
# 从杨图$[\nu]$中按所有可能的方式，去掉一个方块后，所剩下的杨图$[\nu']$，
# 就给出该表示中所包含的$S_{n-1}$群的各种不可约表示（注意$[\nu']$必须也是一个杨图），
# 而且$S_{n-1}$群的每一个不可约表示只出现一次。
#
# 例如：$S_8$群IR[431]中包含$S_7$群IR[43]，[421]，[331]各一次。


# 我们同样关注各种符号、公式、定义背后的物理意义！
# 分支律的意义是，描述可Sn群与Sn-1群的关系
# 不仅仅是+关系，很多具体的量，更是暗藏着带有数学运算的"生成"关系


import copy
import time
from conf.cgc_config import default_s_n, min_s_n_of_branching_law
from core.cgc_utils.cgc_db_typing import BranchingLawInfo
from core.cgc_utils.cgc_local_db import get_branching_laws_finish_s_n_name, get_branching_laws_file_name
from core.young_diagrams import load_young_diagrams
from utils.log import get_logger


logger = get_logger(__name__)


def create_branching_laws(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储分支律实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_branching_law:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_branching_law)
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_branching_laws get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_branching_laws_finish_s_n()
    if not flag:
        err_msg = "get branching_laws finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} branching_laws had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc branching_laws s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的所有杨图的分支律并储存（双循环：外层s_i，内层yd_i）
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        # 这里必须被前面的节点计算过，否则直接推错，而不能再补充计算了
        flag, young_diagrams = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db for calc branching laws meet error with s_i={}, msg={}".format(
                s_i, young_diagrams)
            logger.error(err_msg)
            return False, err_msg
        if not isinstance(young_diagrams, list):
            err_msg = "young_diagrams={} for s_i={} must be list but not".format(young_diagrams, s_i)
            logger.error(err_msg)
            return False, err_msg

        for yd_i in young_diagrams:  # 循环体为[nu_1, nu_2, ...]
            start_time = time.time()
            flag, branching_law_i = calc_single_branching_law(yd_i)  # branching_law_i=(bl_num, rows, cols, before_yd)
            speed_time = int(time.time() - start_time)
            if not flag:
                err_msg = "calc branching_law_i meet error with s_i={}, yd_i={}, msg={}".format(
                    s_i, yd_i, branching_law_i)
                logger.error(err_msg)
                return False, err_msg
            if not branching_law_i or not isinstance(branching_law_i, tuple) or not len(branching_law_i) == 4:
                err_msg = "calc branching_law_i={} should be tuple and len=4".format(branching_law_i)
                logger.error(err_msg)
                return False, err_msg
            bl_num, rows, cols, before_yd = branching_law_i

            # 既然没问题了，那就存（别忘了也要更新Finish_Sn）
            flag, msg = save_single_branching_law(s_i, yd_i, bl_num, rows, cols, before_yd, speed_time)
            if not flag:
                err_msg = "save s_i_branching_law meet error with s_i={}, msg={}".format(s_i, msg)
                logger.error(err_msg)
                return False, err_msg

        s_i_speed_time = int(time.time() - s_i_start_time)
        flag, msg = save_branching_laws_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_branching_laws_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_branching_laws s_n from {} to {} done, return True, finish_s_n={}, using time={}s".format(
        finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_single_branching_law(s_n: int, yd: list,
                              bl_num: int, rows: list, cols: list, before_yd: list, speed_time: int):
    """
    分支律的落盘格式为：
    <top_path>/cgc_results/branching_laws_info/Sn/yd.pkl  ->
    {
    "file_name": "Sn/yd",
    "data": {
            "BL_num": bl_num,    # int，分支数
            "rows": rows,        # list(int)，去掉格子的py行号列表
            "cols": cols,        # list(int)，去掉格子的py列号列表
            "before_YD": [ν’]    # list(list(int))，前置杨盘列表
            }
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    yd表示杨图;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1].pkl: {"BL_num": 2, "rows": [1, 0], "cols": [0, 1], "before_YD": [[2], [1, 1]]}
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_branching_law:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_branching_law)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(bl_num, int) or bl_num <= 0:
        err_msg = "bl_num={} with type={} must be int and > 0".format(bl_num, type(bl_num))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(rows, list):
        err_msg = "rows={} with type={} must be list".format(rows, type(rows))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(cols, list):
        err_msg = "cols={} with type={} must be list".format(cols, type(cols))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(before_yd, list):
        err_msg = "before_yd={} with type={} must be list".format(before_yd, type(before_yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = BranchingLawInfo(s_n)
    _, file_name = get_branching_laws_file_name(s_n, yd)
    table = {"file_name": file_name,
             "data": {
                 "BL_num": bl_num,  # int，分支数
                 "rows": rows,  # list(int)，去掉格子的py行号列表
                 "cols": cols,  # list(int)，去掉格子的py列号列表
                 "before_YD": before_yd  # list(list(int))，前置杨盘列表
             },
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table)
    if not flag:
        return flag, msg

    return True, None


def save_branching_laws_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_branching_law:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_branching_law)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_branching_laws_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = BranchingLawInfo(0)
    _, finish_file_name = get_branching_laws_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {
                           "S{}".format(s_n): s_n_speed_time
                       }}}
    if finish_s_n_before == min_s_n_of_branching_law - 1:
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


def get_branching_laws_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n=0
    """
    _, finish_file_name = get_branching_laws_finish_s_n_name()
    flag, data = BranchingLawInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get branching_laws finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_branching_law - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_branching_law:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_branching_law, data)
        return False, err_msg


def load_branching_law(s_n: int, yd: list, is_flag_true_if_not_s_n=True, is_return_tuple=False):
    """
    取得s_n下young_diagram的分支律(len4tuple，BL_num, rows, cols, before_YD)
    如果没有，根据is_flag_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_branching_law:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_branching_law)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_branching_laws_file_name(s_n, yd)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = BranchingLawInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query branching law with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        bl = data.get("data")
        if bl:
            if is_return_tuple:
                return True, (bl.get("BL_num"), bl.get("rows"), bl.get("cols"), bl.get("before_YD"),)  # bingo！
            else:
                return True, bl  # bingo！

        else:
            err_msg = "branching_law queried from db, but cannot get branching_law from data" \
                      "with data={}, branching_law={} from db".format(data, bl)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist branching_law db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def calc_single_branching_law(young_diagram: list):
    """
    计算分支律
    输入：
    young_diagram====>表示一个杨图ri(list(int))        # [2, 2]
    输出：
    bl_num===========>表示分支个数(int)                # 1
    rows=============>方格所在行号的列表(list(int))     # [1]
    cols=============>方格所在列号的列表(list(int))     # [1]
    before_yd========>前置杨盘(list(list(int)))       # [[2, 1]]
    返回：
    flag, data/msg
    注意：
    这个函数不仅能计算正规杨图，也能计算非正规杨图（如[1, 5, 2]）
    """
    if not isinstance(young_diagram, list):  # 只做粗检查，暂时不用any判断里面必须int了
        err_msg = "young_diagram={} with type={} must be list".format(young_diagram, type(young_diagram))
        logger.error(err_msg)
        return False, err_msg

    if young_diagram == [1]:
        return True, (1, [0], [0], [[]],)

    try:
        yd = copy.deepcopy(young_diagram)
        bl_num = 1
        rows = []
        cols = []
        before_yd = []
        length = len(yd)
        if length < 1:
            err_msg = "len={} of yd={} must >= 1".format(length, yd)
            logger.error(err_msg)
            return False, err_msg
    
        # 计算倒数第一行的分支律
        last_row = length - 1
        if not isinstance(yd[last_row], int) or yd[last_row] < 1:
            err_msg = "yd[{}]={} should be int and >= 1, with yd={}".format(last_row, yd[last_row], yd)
            logger.error(err_msg)
            return False, err_msg
        if yd[last_row] == 1:
            rows += [last_row]
            cols += [0]
            before_yd += [yd[:-1]]
        else:
            yd[last_row] += -1
            rows += [last_row]
            cols += [yd[last_row]]
            before_yd += [yd[:]]
            yd[last_row] += 1

        # 计算倒数第二行到正数第一行的分支律
        if length >= 2:  # 有倒数第二行才有计算
            for i in range(length - 2, -1, -1):  # 从倒数第二行直到第一行 # [l-2, l-3, ..., 0]
                if not isinstance(yd[i], int) or yd[i] < 1:
                    err_msg = "yd[{}]={} should be int and >= 1, with yd={}".format(i, yd[i], yd)
                    logger.error(err_msg)
                    return False, err_msg
                if yd[i] > yd[i + 1]:  # 本行比下面一行多，才能合法去掉一个格子
                    yd[i] += -1  # 真去掉一个，方便后面操作
                    rows += [i]
                    cols += [yd[i]]
                    before_yd += [yd[:]]
                    bl_num += 1  # 计数器加一
                    yd[i] += 1  # 还原成原始的r，再开始下次循环
    except Exception as e:
        err_msg = "calc Branching Law with Young Diagram={} meet error".format(young_diagram)
        logger.error(err_msg)
        logger.error(e)
        return False, e

    # last check
    if not bl_num == len(rows) == len(cols) == len(before_yd):
        err_msg = "this 4 params should be eq but not, " \
                  "its bl_num={} == len(rows={}) == len(cols={}) == len(before_yd={})".format(
            bl_num, rows, cols, before_yd)
        logger.error(err_msg)
        return False, err_msg

    return True, (bl_num, rows, cols, before_yd)
