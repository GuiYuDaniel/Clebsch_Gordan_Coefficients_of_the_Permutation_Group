# -*- coding: utf-8 -*-
"""
this code for creating Young Tableaux by Sn, young_diagrams and branching_laws
计算Sn杨图的杨盘
"""

# \subsection{杨盘}
#
# 在$S_n$群杨图$[\nu]$上按下面规则填上1，2，...，n 个数字所构成的一个图称为杨盘。
# 规则要求我们在任一行自左向右和任一列自上向下阅读时，这些数字必须是按增加次序排列的。
# 我们用$Y^{[\nu]}_m$来代表一个杨盘，m为杨盘的编号。
#
# 我们可用杨盘来标志$S_n$的标准基$|Y^{[\nu]}_m>$。
# 它的意义是按照这样来规定的：
# 在杨盘$Y^{[\nu]}_m$中去掉方块n，就得到涉及排列好的n-1个方块的一个盘$Y^{[\nu']}_{m'}$；
# 再去掉方块n-1，就又得到排列好的下一个盘$Y^{[\nu'']}_{m''}$，等等，直到最后。
# 规定不可约的基失$|Y^{[\nu]}_m>$是这样一种基失：
# 它分别属于群$S_n,S_{n-1},S_{n-2},...,S_1$的IR$[\nu],[\nu'],[\nu''],...,[1]$。


# 我们同样关注各种符号、公式、定义背后的物理意义！
# 利用杨图/杨盘可以直观地分类置换群格子之间的对称性（详见《群表示论的新途径》陈金全（上海科学技术出版社1984）2.6；3.4）
# 比如我们看到 [[1, 2, 3]]这个杨盘，
#   a）当它表示粒子的时候，可以代表1，2，3之间具有交换不变性，比方说全同粒子
#   b）当它表示波函数的时候，
#       b1）可以表示α^3（每个粒子可能取值只有α这1种状态）既三个全同粒子，的全部对称性状态(|ααα>)
#       b2）可以表示(α^2)β（每个粒子可能取值有α和β这2种状态，但有两个是一致的）的部分对称状态（群空间，不是坐标空间，是表示三个群空间基失全对称）
#           （说明：坐标空间可以用|ααβ>，|αβα>，|βαα>表示，但群空间上他们体现不出对称性关系；
#                 群空间里重新线性组合并正交归一化，
#                 用(|ααβ>+|αβα>+|βαα>)/sqrt(3), (2|ααβ>-|αβα>-|βαα>)/sqrt(6), (|αβα>-|βαα>)/sqrt(2) 三个新基失表示，
#                 其中[[1, 2, 3]]就代表(|ααβ>+|αβα>+|βαα>)/sqrt(3)。
#                 这种变换后不仅具有对称性清晰一个好处，而且本征值也便于计算）
#       b3）可以表示αβγ（每个粒子可能取值有α、β、γ这3种状态，可能一致也可能不一致，都存在）的六种状态中的，全同那个状态


import copy
import time
from core.branching_laws import load_branching_law
from conf.cgc_config import default_s_n, min_s_n_of_young_table
from core.cgc_utils.cgc_local_db import get_young_tableaux_file_name, get_young_tableaux_finish_s_n_name, \
    get_young_tableaux_num_file_name, get_young_tableaux_phase_factor_file_name
from core.cgc_utils.cgc_db_typing import YoungTableInfo
from core.young_diagrams import load_young_diagrams
from core.young_diagrams import is_young_diagram, calc_s_n_from_young_diagram, calc_s_n_from_young_diagram_without_check
from functools import lru_cache
from utils.log import get_logger


logger = get_logger(__name__)


def create_young_tableaux(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储杨盘(以及杨盘总数)实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_young_tableaux get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_young_tableaux_finish_s_n()
    if not flag:
        err_msg = "get young_tableaux finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} young_tableaux had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc young_tableaux s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的所有杨图的杨盘并储存（双循环：外层s_i，内层yd_i）
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        s_i_start_time = time.time()
        # 这里必须被前面的节点计算过，否则直接推错，而不能再补充计算了
        flag, young_diagrams = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)  # True是给cmd的
        if not flag:
            err_msg = "get young_diagrams db for calc young_tableaux meet error with s_i={}, msg={}".format(
                s_i, young_diagrams)
            logger.error(err_msg)
            return False, err_msg
        if not isinstance(young_diagrams, list):
            err_msg = "young_diagrams={} for s_i={} must be list but not".format(young_diagrams, s_i)
            logger.error(err_msg)
            return False, err_msg

        for yd_i in young_diagrams:  # 循环体为[nu_1, nu_2, ...]
            start_time = time.time()
            # 计算young_table
            flag, young_table_i = calc_single_young_table(s_i, yd_i)  # young_table_i={"m_j": young_table_j}
            if not flag:
                err_msg = "calc young_table_i meet error with s_i={}, yd_i={}, msg={}".format(s_i, yd_i, young_table_i)
                logger.error(err_msg)
                return False, err_msg
            if not young_table_i or not isinstance(young_table_i, dict):
                err_msg = "calc young_table_i={} should be dict".format(young_table_i)
                logger.error(err_msg)
                return False, err_msg

            # 计算young_table的phase_factor，按照顺序存列表
            flag, phase_factor_list = quick_calc_phase_factor_list(s_i, yd_i, young_table_i)
            if not flag:
                err_msg = "calc phase_factor_list meet error with s_i={}, yd_i={}, msg={}".format(
                    s_i, yd_i, young_table_i)
                logger.error(err_msg)
                return False, err_msg

            speed_time = int(time.time() - start_time)
            # 既然没问题了，那就存（别忘了也要更新Finish_Sn）
            flag, msg = save_single_young_table(s_i, yd_i, young_table_i, speed_time)
            if not flag:
                err_msg = "save young_table_i meet error with s_i={}, yd_i={}, msg={}".format(
                    s_i, yd_i, msg)
                logger.error(err_msg)
                return False, err_msg

            # young_table特别另存一个young_table_num
            flag, msg = save_single_young_table_num(s_i, yd_i, len(young_table_i))
            if not flag:
                err_msg = "save save_single_young_table_num meet error with s_i={}, yd_i={}, msg={}".format(
                    s_i, yd_i, msg)
                logger.error(err_msg)
                return False, err_msg

            # young_table特别另存一个young_table_Λ
            flag, msg = save_single_young_table_phase_factor(s_i, yd_i, phase_factor_list)
            if not flag:
                err_msg = "save save_single_young_table_phase_factor meet error with s_i={}, yd_i={}, msg={}".format(
                    s_i, yd_i, msg)
                logger.error(err_msg)
                return False, err_msg

        # 别忘了也要更新Finish_Sn
        s_i_speed_time = int(time.time() - s_i_start_time)
        flag, msg = save_young_tableaux_finish_s_n(s_i, s_i_speed_time, is_check_add_one=True)
        if not flag:
            err_msg = "save_young_tableaux_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    c_time = time.time() - start_time_c
    logger.info("#### create_young_tableaux s_n from {} to {} done, return True, finish_s_n={}, using time={}s".format(
        finish_s_n + 1, s_n, s_n, c_time))
    return True, s_n


def save_single_young_table(s_n: int, yd: list, young_table: dict, speed_time: int):
    """
    杨盘的落盘格式为：
    <CG>/young_tableaux_info/Sn/[ν_i].pkl  ->
    {
    "file_name": "Sn/[ν_i]",
    "data": {"m_i": young_table}
    "flags": {"speed_time": speed_time
              "total_num": total_num}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    total_num就是len({"m_i": young_table})，表示杨盘总数;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1].pkl: {"1": [[1, 2], [3]], "2": [[1, 3], [2]]}
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(young_table, dict):
        err_msg = "yd={} with type={} must be dict".format(young_table, type(young_table))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg

    db_info = YoungTableInfo(s_n)
    _, file_name = get_young_tableaux_file_name(s_n, yd)
    total_num = len(young_table)
    table = {"file_name": file_name,
             "data": young_table,
             "flags": {"speed_time": speed_time,
                       "total_num": total_num}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table, point_key="data")
    if not flag:
        return flag, msg

    return True, None


def save_single_young_table_num(s_n: int, yd: list, young_table_num: int):
    """
    杨盘总数的落盘格式为：
    <CG>/young_tableaux_info/Sn/[ν_i]_num.pkl  ->
    {
    "file_name": "Sn/[ν_i]_num",
    "data": {"total_num": total_num}
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    total_num就是len({"m_i": young_table})，表示杨盘总数;

    例如：
    S3/[2, 1]_num.pkl: {"total_num": 2}
    """
    # 省略检查 前面都查过了

    db_info = YoungTableInfo(s_n)
    _, file_name = get_young_tableaux_num_file_name(s_n, yd)
    table = {"file_name": file_name,
             "data": young_table_num,
             "flags": {}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table, point_key="data")
    if not flag:
        return flag, msg

    return True, None


def save_single_young_table_phase_factor(s_n: int, yd: list, young_table_phase_factor_list: list):
    """
    杨盘相位因子的落盘格式为：
    <CG>/young_tableaux_info/Sn/[ν_i]_Λ.pkl  ->
    {
    "file_name": "Sn/[ν_i]_Λ",
    "data": phase_factor_list
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    phase_factor_list就是Λ按照m从小到大，或者说Yamanouchi序排列的

    例如：
    S3/[2, 1]_Λ.pkl: [1, -1]
    """
    # 检查
    if not isinstance(young_table_phase_factor_list, list):
        err_msg = "young_table_phase_factor_list={} with type={} must be list".format(
            young_table_phase_factor_list, type(young_table_phase_factor_list))
        logger.error(err_msg)
        return False, err_msg

    db_info = YoungTableInfo(s_n)
    _, file_name = get_young_tableaux_phase_factor_file_name(s_n, yd)
    table = {"file_name": file_name,
             "data": young_table_phase_factor_list,
             "flags": {}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table, point_key="data")
    if not flag:
        return flag, msg

    return True, None


def save_young_tableaux_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_young_tableaux_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = YoungTableInfo(0)
    _, finish_file_name = get_young_tableaux_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": {},
             "flags": {"finish_s_n": s_n,
                       "history_times": {
                           "S{}".format(s_n): s_n_speed_time
                       }}}
    if finish_s_n_before == min_s_n_of_young_table - 1:
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


def get_young_tableaux_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n=0
    """
    _, finish_file_name = get_young_tableaux_finish_s_n_name()
    flag, data = YoungTableInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get young_tableaux finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, min_s_n_of_young_table - 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= min_s_n_of_young_table:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and >= {}, with data={}".format(finish_s_n, min_s_n_of_young_table, data)
        return False, err_msg


def calc_single_young_table(s_n: int, young_diagram: list):
    """
    计算杨盘
    输入：
    s_n==============>表示当前置换群阶数(int)                   # 3
    young_diagram====>表示一个杨图ri(list(int))                # [2, 1]
    输出：
    young_table======>杨图ri的所有杨盘(dict(list(list(int))))  # {"1":[[1,2],[3]], "2":[[1,3],[2]]}
    返回：
    flag, data/msg
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(young_diagram, list):
        err_msg = "young_diagram={} with type={} must be list".format(young_diagram, type(young_diagram))
        logger.error(err_msg)
        return False, err_msg
    if s_n != sum(young_diagram):
        err_msg = "sum_young_diagram={} must eq s_n={} but not, with young_diagram={}".format(
            sum(young_diagram), s_n, young_diagram)
        logger.error(err_msg)
        return False, err_msg

    # 保护
    if s_n == 1:
        return True, {"1": [[1]]}

    # 取前置数据
    flag, bl_tuple = load_branching_law(s_n, young_diagram, is_flag_true_if_not_s_n=False, is_return_tuple=True)
    if not flag:
        err_msg = "get branching_law db for calc young_table meet error with s_n={}, msg={}".format(s_n, bl_tuple)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(bl_tuple, tuple) or len(bl_tuple) != 4:
        err_msg = "bl_tuple={} for s_n={}, young_diagram={} must be tuple and len 4 but not".format(
            bl_tuple, s_n, young_diagram)
        logger.error(err_msg)
        return False, err_msg
    bl_num, rows, cols, before_yd = bl_tuple
    if not isinstance(bl_num, int) or bl_num <= 0:
        err_msg = "bl_num={} with type={} must be int and > 0".format(bl_num, type(bl_num))
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
    s_b = s_n - 1  # Sn_before 如果多个 就叫 s_b_1, s_b_2, ...

    # 真正开始计算杨盘
    # 双循环：
    # 外层是按照当前young_diagram分支律的分支数bl_num开启循环，其中，range(bl_num)与rows,cols,before_yd一一对应
    # 内层是按照当前分支before_young_table_i_batch开启循环，循环节是s_b的单个杨盘py序，对Sn的意义是block内部的序号
    # 希望i标志对于Sn的、整体的序号；j标志对于Sb的、块内的序号
    total_num = 0
    young_table = {}
    for bl_index, row_i, col_i, before_yd_i in zip(range(bl_num), rows, cols, before_yd):
        flag, before_young_table_i_batch = load_young_table(s_b, before_yd_i)  # 第i个分支的杨盘(Sn-1的)
        if not flag:
            err_msg = "get before_young_table_i_batch by s_b={}, before_yd_i={} meet error with msg={}".format(
                s_b, before_yd_i, before_young_table_i_batch)
            logger.error(err_msg)
            return False, err_msg
        before_young_table_i_batch_total_num = len(before_young_table_i_batch)  # 必须先手读取杨盘的时候，len计算杨盘总数就可以了

        # 这是分支律中的第一个，一定是倒数第一行的格子(不一定是1哦)被去掉。分两种情况：
        # 分支律中的第一个，且去掉位是倒数第一行、正数第一列（特殊性是填在Sn格子的时候要新起一行）
        if bl_index == 0 and col_i == 0:
            for block_j in range(before_young_table_i_batch_total_num):  # 对每个before构型都要分别建立新的杨盘
                math_block_j = block_j + 1
                math_total_i = math_block_j
                before_yt_j = before_young_table_i_batch.get(str(math_block_j))  # Sn-1的第j个杨盘
                before_yt_j.append([s_n])  # 生成Sn杨图，通过在尾部添加一行且一个格子
                young_table[str(math_total_i)] = before_yt_j

        # 分支律的非第一分支，或去掉位不是倒数第一行、正数第一列（无需添加新行，但计算math_block_j时要考虑到前面分支的贡献）
        else:
            for block_j in range(before_young_table_i_batch_total_num):
                math_block_j = block_j + 1
                math_total_i = total_num + math_block_j
                before_yt_j = before_young_table_i_batch.get(str(math_block_j))  # Sn-1的第j个杨盘
                before_yt_j[row_i].append(s_n)  # 生成Sn杨图，通过在row_i行尾添加一个格子
                young_table[str(math_total_i)] = before_yt_j

        total_num += before_young_table_i_batch_total_num

    # double check
    if len(young_table) != total_num:
        err_msg = "young_table_len={} must eq total_num={}, pls check, with s_n={}, young_diagram={}".format(
            len(young_table), total_num, s_n, young_diagram)
        logger.error(err_msg)
        return False, err_msg

    return True, young_table


def load_young_table(s_n: int, yd: list, is_flag_true_if_not_s_n=True):
    """
    取得s_n下young_diagram的杨盘(字典：键是杨盘编号(str)；值是编号对应的杨盘(list(list(int))))
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_young_tableaux_file_name(s_n, yd)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = YoungTableInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query young table with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        young_table = data.get("data")
        if young_table and isinstance(young_table, dict):
            # 只检查有没有 不对内容做检查了
            return True, young_table  # bingo！

        else:
            err_msg = "young_table queried from db, but cannot get young_table from data" \
                      "with data={}, young_table={} from db".format(data, young_table)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist young_table db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def load_young_table_num(s_n: int, yd: list, is_flag_true_if_not_s_n=True):
    """
    取得s_n下young_diagram的杨盘总数
    如果没有，根据is_flag_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_young_tableaux_num_file_name(s_n, yd)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = YoungTableInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query young_table_num with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        young_table_num = data.get("data", None)
        if young_table_num and isinstance(young_table_num, int):
            # 只检查有没有 不对内容做检查了
            return True, young_table_num  # bingo！

        else:
            err_msg = "data queried from db, but cannot get young_table_num from data" \
                      "with data={}, young_table_num={} from db".format(data, young_table_num)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist young_table_num db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def load_young_table_phase_factor(s_n: int, yd: list, is_flag_true_if_not_s_n=True):
    """
    取得s_n下young_diagram的相位因子列表
    如果没有，根据is_flag_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_young_tableaux_phase_factor_file_name(s_n, yd)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = YoungTableInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query young_table_num with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        young_table_phase_factor = data.get("data", None)
        if young_table_phase_factor and isinstance(young_table_phase_factor, list):
            # 只检查有没有 不对内容做检查了
            return True, young_table_phase_factor  # bingo！

        else:
            err_msg = "data queried from db, but cannot get young_table_phase_factor from data" \
                      "with data={}, young_table_phase_factor={} from db".format(data, young_table_phase_factor)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist young_table_phase_factor db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


def is_young_table(young_table):
    """
    1，检查类型为list
    2，检查非空
    3，检查内部类型为list非空
    4，检查内部的内部为int
    5，检查升序性
    6，检查满编
    """
    if not young_table or not isinstance(young_table, list):
        err_msg = "young_table={} with type={} must be real list".format(young_table, type(young_table))
        logger.error(err_msg)
        return False

    yd = []
    all_number = []
    for yt_layer in young_table:
        if not yt_layer or not isinstance(yt_layer, list):
            err_msg = "yt_layer={} with type={} must be real list".format(yt_layer, type(yt_layer))
            logger.error(err_msg)
            return False
        if not all(isinstance(i, int) for i in yt_layer):
            err_msg = "yt_layer={} must have all int element but {}".format(yt_layer, [type(i) for i in yt_layer])
            logger.error(err_msg)
            return False
        ascending_yt_layer = sorted(yt_layer)
        if yt_layer != ascending_yt_layer:
            err_msg = "yt_layer={} must be ascending order, like {}".format(yt_layer, ascending_yt_layer)
            logger.error(err_msg)
            return False
        yd.append(len(yt_layer))
        all_number += yt_layer

    if not is_young_diagram(yd):
        err_msg = "young_table={} have wrong young_diagram={}, pls check".format(young_table, yd)
        logger.error(err_msg)
        return False

    ascending_all_number = sorted(all_number)
    if ascending_all_number != list(range(1, len(all_number) + 1)):
        err_msg = "young_table={} have wrong all_number={}, it must be {}".format(
            young_table, all_number, list(range(1, len(all_number) + 1)))
        logger.error(err_msg)
        return False

    return True


def calc_young_diagram_from_young_table_without_check(young_table):
    yd = [len(i) for i in young_table]
    return True, yd


def calc_young_diagram_from_young_table(young_table):
    if not is_young_table(young_table):
        err_msg = "cannot calc young_diagram with wrong young_table={}".format(young_table)
        logger.error(err_msg)
        return False, None
    return calc_young_diagram_from_young_table_without_check(young_table)


def calc_s_n_from_young_table_without_check(young_table):
    s_n = max([max(i) for i in young_table])
    return True, s_n


def calc_s_n_from_young_table(young_table):
    if not is_young_table(young_table):
        err_msg = "cannot calc s_n with wrong young_table={}".format(young_table)
        logger.error(err_msg)
        return False, None

    return calc_s_n_from_young_table_without_check(young_table)


def get_s_i_index_in_young_table(s_i, young_table):
    """返回s_i格子在杨盘中的py坐标"""
    yt_row_num = len(young_table)
    for row, sub_list in zip(range(yt_row_num), young_table):
        if s_i in sub_list:
            col = sub_list.index(s_i)
            return row, col

    return False, "s_i={} not in young_table={}, pls check".format(s_i, young_table)


def _quickly_calc_in_decreasing_page_order():
    """
    TODO 递归计算，可以依赖已有结果，也可以不依赖 如果需要 把它完成并且和不带_的函数交换
    """
    pass


def read_young_table_in_decreasing_page_order(young_table, all_yts_with_same_yt, is_rst_int=False):
    """
    1，错误：返回False, err_msg
    2，没有：返回True, None
    3，有  ：返回True, str/int(idp_order)
    """
    # 简单检查
    if not young_table or not isinstance(young_table, list):
        err_msg = "young_table={} with type={} must be real list".format(young_table, type(young_table))
        logger.error(err_msg)
        return False, err_msg
    if not all_yts_with_same_yt or not isinstance(all_yts_with_same_yt, dict):
        err_msg = "all_yts_with_same_yt={} with type={} must be real dict".format(
            all_yts_with_same_yt, type(all_yts_with_same_yt))
        logger.error(err_msg)
        return False, err_msg

    # 判断
    for idp_order, v in all_yts_with_same_yt.items():
        if young_table == v:
            if is_rst_int:
                return True, int(idp_order)
            else:
                return True, idp_order
    return True, None


def quickly_calc_young_table_in_decreasing_page_order(young_table, yd=None, s_n=None, is_rst_int=False):
    """
    快速计算杨盘编号（也称作Yamanouchi编号）（不递归）
    TODO 这就是我那个定理嘛 1，独立成可以对外的api；2，贴上描述
    """
    # 检查手续
    if not is_young_table(young_table):
        err_msg = "young_table={} with must be is_young_table".format(young_table)
        logger.error(err_msg)
        return False, err_msg
    if yd is None:
        flag, yd = calc_young_diagram_from_young_table_without_check(young_table)
        if not flag:
            err_msg = "calc_young_diagram_from_young_table_without_check meet error with msg={}".format(yd)
            logger.error(err_msg)
            return False, err_msg
    else:
        # if not is_young_diagram(yd):
        #     err_msg = "yd={} with must be is_young_diagram".format(yd)
        #     logger.error(err_msg)
        #     return False, err_msg
        if not yd or not isinstance(yd, list):  # 对于传入的yd为了快，我们只简单检查
            err_msg = "yd={} with type={} must be real list".format(yd, type(yd))
            logger.error(err_msg)
            return False, err_msg
    if s_n is None:
        flag, s_n_from_yt = calc_s_n_from_young_table_without_check(young_table)
        if not flag:
            err_msg = "calc_s_n_from_young_table_without_check meet error with msg={}".format(s_n_from_yt)
            logger.error(err_msg)
            return False, err_msg
        flag, s_n_from_yd = calc_s_n_from_young_diagram_without_check(yd)
        if not flag:
            err_msg = "calc_s_n_from_young_diagram_without_check meet error with msg={}".format(s_n_from_yd)
            logger.error(err_msg)
            return False, err_msg
        if s_n_from_yt != s_n_from_yd:
            err_msg = "calc s_n={} by yt={} must eq s_n={} by yd={}".format(s_n_from_yt, young_table, s_n_from_yd, yd)
            logger.error(err_msg)
            return False, err_msg
        else:
            s_n = s_n_from_yt
    else:
        if not s_n >= 1 or not isinstance(s_n, int):
            err_msg = "s_n={} with type={} must be >=1 and int".format(s_n, type(s_n))
            logger.error(err_msg)
            return False, err_msg

    # 计算程序
    idp_order = 0
    remain_yd_s_max = copy.deepcopy(yd)
    for s_max in range(s_n, 0, -1):  # [Sn, Sn-1, Sn-2, ..., 1]  # s_max表示剩余杨图的最大格子
        if s_max == 1:
            idp_order += 1
        else:
            flag, bl_tuple_s_max = load_branching_law(s_max, remain_yd_s_max, is_flag_true_if_not_s_n=False,
                                                      is_return_tuple=True)
            if not flag:
                err_msg = "get branching_law for calc quickly_calc_young_table_in_decreasing_page_order meet error " \
                          "with s_max={}, remain_yd_s_max={}, msg={}".format(s_max, remain_yd_s_max, bl_tuple_s_max)
                logger.error(err_msg)
                return False, err_msg
            bl_num_s_max, _, _, before_yd_s_max = bl_tuple_s_max

            # 计算去掉当前最大格子的剩余杨图
            row_s_max, _ = get_s_i_index_in_young_table(s_max, young_table)
            if remain_yd_s_max[row_s_max] == 1:  # 表示剩余最大格子在剩余杨图中那行只剩一个格子了，该行就一定是尾行
                remain_yd_s_max = remain_yd_s_max[:-1]  # 去掉这种格子，相当于去掉尾行
            else:  # 表示剩余最大格子在剩余杨图中那行不止一个格子
                remain_yd_s_max[row_s_max] += -1  # 去掉这种格子，相当于对该行减1

            # 一定要在loadBeforeR之后再更新删去sn的杨盘
            for bl_j_s_max, bl_j_s_max_yd in zip(range(bl_num_s_max), before_yd_s_max):  # 对前面分支能容纳的杨盘总数求和
                if bl_j_s_max_yd == remain_yd_s_max:
                    # 终止条件，它就是去掉s_max后的yd那个分支，从包括这个分支开始以及后面的分支，对杨盘编号没有贡献
                    break
                flag, before_yt_s_max_j_num = load_young_table_num(s_max-1, bl_j_s_max_yd,
                                                                   is_flag_true_if_not_s_n=False)
                if not flag:  # 这里简化检查手续了
                    err_msg = "get before_yt_s_max_j_num for calc quickly_calc_young_table_in_decreasing_page_order " \
                              "by s_max-1={}, bl_j_s_max_yd={} meet error " \
                              "with msg={}".format(s_max-1, bl_j_s_max_yd, before_yt_s_max_j_num)
                    logger.error(err_msg)
                    return False, err_msg
                idp_order += before_yt_s_max_j_num

    if is_rst_int:
        return True, idp_order
    else:
        return True, str(idp_order)


def _calc_single_phase_factor(y_target, y_start, is_output_exchange_list=False):
    """
    使用下面的算法（优化），总能把任意y_start通过有限次临近交换，变换成y_target：
    # 很像最笨的方法玩魔方，每次轮只把一个目标数字移动到最终位置，下一轮就可以不触动它和更小的数字，去移动更大的数字了
    1，找不同（y_xxx 不使用标准yt，使用list(int)就够了）
    2，得到不同中最小的数
    3，该数字出现在目标中的位置
    4，该位置，对应start中的数字
    5，使用临近交换，把y_start中上述最小数字换过去（不用逐个交换，直接使用(k-1, k)(k-2, k-1)...(i,i+1)的结果就可以了）
    （结果就是：i去到k位置，其他数字，顺序去到前一个数字的位置；使用了k-i次交换）
    6，得到交换后的y_start
    7，循环，直到y_target == y_start

    对于is_output_exchange_list开关，常规是关闭的，但如果有需要，可以通过传入True，得到形如[(4,5),(3,4),(2,3),(5,6)]这样的具体交换细节
    注意：is_output_exchange_list会改变输出个数

    除了直接使用定义，还可以用Butler的方法，用(-1)^nb，乘以对应分支的Λ_st，得到Λ
    """
    exchange_times = 0
    exchange_list = []  # 具体的交换顺序
    while y_target != y_start:
        diff_num_list = [start for target, start in zip(y_target, y_start) if target != start]
        min_diff_num = min(diff_num_list)  # 由于diff_num_list生成方式，也是diff_num_list[0]
        target_index_of_min_diff_num = y_target.index(min_diff_num)
        start_num_of_target_index = y_start[target_index_of_min_diff_num]
        # 开始通过改造start生成end（注意，由于python机制，y_start会变）
        y_start[target_index_of_min_diff_num] = min_diff_num  # 首先把min_diff_num换过去
        # 然后平移剩下的
        for num in range(start_num_of_target_index-1, min_diff_num-1, -1):  # k-1 ~ i
            # 点到名字的数字加一，等价于i~k-1左移动
            # 其中，列表中有两个min_diff_num，只拿第一个index是合理利用python机制
            y_start[y_start.index(num)] += 1
            if is_output_exchange_list:
                exchange_list.append(tuple([num, num+1]))
        exchange_times += start_num_of_target_index - min_diff_num
    if exchange_times % 2 == 0:
        phase_factor = 1  # 按照定义，偶数次交换，相位因子为1
    else:
        phase_factor = -1  # 按照定义，奇数次交换，相位因子为-1

    if is_output_exchange_list:
        return phase_factor, exchange_list
    else:
        return phase_factor


def calc_single_yt_phase_factor(young_table_target, _is_check_yt=True, is_output_exchange_list=False):
    """根据定义的原始方法，遍历临近交换树，计算Λ
    作为独立对外的API，建议检查yt

    对于_is_check_yt开关，默认打开，如需关闭，请保证传入的young_table_target一定正确
    对于is_output_exchange_list开关，默认关闭，但如果有需要，可以通过传入True，得到具体交换细节
    注意：is_output_exchange_list会改变输出个数
    """
    if _is_check_yt and not is_young_table(young_table_target):
        err_msg = "young_table_target={} should be a real young table".format(young_table_target)
        logger.error(err_msg)
        return False, err_msg

    y_target = []
    for yt_row in young_table_target:
        y_target += yt_row
    y_first = list(range(1, len(y_target) + 1))
    # _calc_single_phase_factor不改变young_table_target，只改变y_first，所以是安全的
    return _calc_single_phase_factor(y_target, y_first, is_output_exchange_list=is_output_exchange_list)


def quick_calc_phase_factor_list(s_n: int, yd: list, young_tableaux: dict):
    """使用论文中的推论1，快速计算全部yt的相位因子Λ"""
    if not isinstance(s_n, int) or s_n < min_s_n_of_young_table:
        err_msg = "s_n={} with type={} must be int and >= {}".format(s_n, type(s_n), min_s_n_of_young_table)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "young_diagram={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(young_tableaux, dict):
        err_msg = "young_tableaux={} with type={} must be dict".format(young_tableaux, type(young_tableaux))
        logger.error(err_msg)
        return False, err_msg

    # 保护
    if s_n == 1:
        return True, [1]

    # 取前置数据
    _, bl_tuple = load_branching_law(s_n, yd, is_flag_true_if_not_s_n=False, is_return_tuple=True)
    _, _, _, before_yd = bl_tuple
    s_b = s_n - 1

    # 开始计算
    phase_factor_list = []
    for yd_sb in before_yd:
        leader_yt = young_tableaux.get(str(len(phase_factor_list) + 1))
        leader_yt_phase_factor = calc_single_yt_phase_factor(leader_yt, _is_check_yt=False)
        flag, sb_phase_factor_list = load_young_table_phase_factor(s_b, yd_sb, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "get sb_phase_factor_list meet error with s_i={}, yd={}, msg={}".format(
                s_b, yd_sb, sb_phase_factor_list)
            logger.error(err_msg)
            return False, err_msg
        if leader_yt_phase_factor != sb_phase_factor_list[0]:
            sb_phase_factor_list = [-1 * i for i in sb_phase_factor_list]
        phase_factor_list += sb_phase_factor_list

    return True, phase_factor_list
