# -*- coding: utf-8 -*-
"""
this code for creating Yamanouchi Matrix (ij)(in) by Sn, young_diagrams, branching_laws and young_tableaux
计算Sn杨图的Yamanouchi对换矩阵
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）第四章第4节：置换群标准矩阵元
# Yamanouchi对换矩阵，是一种将群元的IR看为基失和算符，算符作用在基失上，可以将标准基|Yi>线性变换为|Yj>+|Yk>的矩阵
# （或称符合Yamanouchi相位约定的一种IR矩阵）
# 注意区分杨盘Ym和标准基|Ym>
# 当(i-1,i)Ym=Yn时，(i-1,i)|Ym> != |Yn>
# 而是(i-1,i)|Ym> = 1/σ |Ym> + sqrt(σ^2 - 1)/abs(σ) |Yn>
# 其中，σ表示轴距
# 式中，取相邻对换的非对角矩阵元永远大于等于零。这种相位约定称为Young-Yamanouchi标准相位。


# 公式（Hamermesh）
# a）当i-1和i在Sn群的杨盘Y^{[\nu]}_m中处于同一行或同一列时
#    (i-1, i)|Y^{[\nu]}_m> = +- |Y^{[\nu]}_m>    # 注意，两遍都是Y^{[\nu]}_m，并没有改换成非法标准基哦！
# b）当i-1和i不在同一行或列时
#    D^{[\nu]}_nm(i-1, i) = <Y^{[\nu]}_n>|(i-1, i)|Y^{[\nu]}_m>
#                            /   1/σ                    当n=m，
#                         = {   sqrt(σ^2 - 1)/abs(σ)，  当Y^{[\nu]}_n = (i-1, i)Y^{[\nu]}_m 时，（也得是合法杨盘）
#                            \   0                      其他情形
# 矩阵里，就没有非法标准基或者非法杨盘的位置！


# 我们同样关注各种符号、公式、定义背后的物理意义！
# Yamanouchi对换矩阵，又叫置换群标准矩阵元
# 但要注意，通过这个矩阵（二循环/对换），是将一个杨盘转换成新的杨盘，却不是把一个标准基转化成另外一个，因为标准基不直接是杨盘


import copy
import numpy as np
import math
import time
from core.branching_laws import load_branching_law
from conf.cgc_config import default_s_n
from core.cgc_utils.cgc_local_db import get_yamanouchi_matrix_file_name, get_yamanouchi_matrix_finish_s_n_name
from core.cgc_utils.cgc_db_typing import YamanouchiMatrixInfo
from core.young_diagrams import load_young_diagrams
from core.young_tableaux import load_young_table, load_young_table_num
from core.young_tableaux import get_s_i_index_in_young_table
from core.young_tableaux import read_young_table_in_decreasing_page_order
from core.young_tableaux import quickly_calc_young_table_in_decreasing_page_order
from utils.log import get_logger


logger = get_logger(__name__)


def create_yamanouchi_matrix(s_n: int=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储Yamanouchi矩阵实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg

    logger.info("#### create_yamanouchi_matrix get input s_n={}".format(s_n))
    start_time_c = time.time()

    # 先查询数据库中完成到S几：如果输入s_n未计算，直接从循环中cut掉算好的部分；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_yamanouchi_matrix_finish_s_n()
    if not flag:
        err_msg = "get yamanouchi_matrix finish_s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n == 1:
        msg = "s_n={} have no yamanouchi_matrix, jump it and return True"
        logger.info(msg)
        return True, s_n
    elif s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} yamanouchi_matrix had been calculated, return True, s_n".format(s_n)
        logger.info(msg)
        return True, s_n
    else:
        msg = "finish_s_n={}, will calc yamanouchi_matrix s_n from {} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的所有杨图的杨盘并储存（双循环：外层s_i，内层yd_i）
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        start_time_s_i = time.time()
        # 这里必须被前面的节点计算过，否则直接推错，而不能再补充计算了
        flag, young_diagrams = load_young_diagrams(s_i, is_flag_true_if_not_s_n=False)  # True是给cmd的
        if not flag:
            err_msg = "get young_diagrams db for calc yamanouchi_matrix meet error with s_i={}, msg={}".format(
                s_i, young_diagrams)
            logger.error(err_msg)
            return False, err_msg
        if not isinstance(young_diagrams, list):
            err_msg = "young_diagrams={} for s_i={} must be list but not".format(young_diagrams, s_i)
            logger.error(err_msg)
            return False, err_msg

        for yd_i in young_diagrams:  # 循环体为[nu_1, nu_2, ...]
            start_time_3 = time.time()

            key_ij_list = [(s_i - i - 1, s_i - i,) for i in range(s_i - 1)]
            key_in_list = [(s_i - i - 1, s_i,) for i in range(s_i - 1)]
            for key_ij, key_in, rst_iter in zip(key_ij_list, key_in_list, calc_yamanouchi_matrix_iter(s_i, yd_i)):
                flag, ym_ij_dict, ym_in_dict = rst_iter
                if not flag:
                    err_msg = "calc yamanouchi_matrix_ix meet error with s_i={}, yd_i={}, msg={}".format(
                        s_i, yd_i, ym_ij_dict)
                    logger.error(err_msg)
                    return False, err_msg
                if not ym_ij_dict or not isinstance(ym_ij_dict, dict) \
                        or not isinstance(ym_ij_dict.get(key_ij), np.ndarray):
                    err_msg = "calc ym_ij_dict={} meet error with key_ij={}, yd_i={}".format(ym_ij_dict, key_ij, yd_i)
                    logger.error(err_msg)
                    return False, err_msg
                if not ym_in_dict or not isinstance(ym_in_dict, dict) \
                        or not isinstance(ym_in_dict.get(key_in), np.ndarray):
                    err_msg = "calc ym_in_dict={} meet error with key_in={}, yd_i={}".format(ym_in_dict, key_in, yd_i)
                    logger.error(err_msg)
                    return False, err_msg
                ym_ij = ym_ij_dict.get(key_ij)
                ym_in = ym_in_dict.get(key_in)

                # 既然没问题了，那就存
                speed_time_3 = int(time.time() - start_time_3)
                flag, msg = save_yamanouchi_matrix(s_i, yd_i, key_ij, ym_ij, speed_time_3, mode="ij")
                if not flag:
                    err_msg = "save yamanouchi_matrix_ij meet error " \
                              "with s_i={}, yd_i={}, key_ij={}, ym_ij={}, speed_time_3={}, " \
                              "msg={}".format(s_i, yd_i, key_ij, ym_ij, speed_time_3, msg)
                    logger.error(err_msg)
                    return False, err_msg
                flag, msg = save_yamanouchi_matrix(s_i, yd_i, key_in, ym_in, speed_time_3, mode="in")
                if not flag:
                    err_msg = "save yamanouchi_matrix_in meet error " \
                              "with s_i={}, yd_i={}, key_in={}, ym_in={}, speed_time_3={}, " \
                              "msg={}".format(s_i, yd_i, key_in, ym_in, speed_time_3, msg)
                    logger.error(err_msg)
                    return False, err_msg
                start_time_3 = time.time()

        # 别忘了也要更新Finish_Sn
        speed_time_s_i = int(time.time() - start_time_s_i)
        flag, msg = save_yamanouchi_matrix_finish_s_n(s_i, speed_time_s_i, is_check_add_one=True)
        if not flag:
            err_msg = "save_yamanouchi_matrix_finish_s_n meet error with s_i={}, msg={}".format(s_i, msg)
            logger.error(err_msg)
            return False, err_msg

    speed_time_c = time.time() - start_time_c
    logger.info("#### create_yamanouchi_matrix s_n from {} to {} done, return True, finish_s_n={}, "
                "using time={}s".format(finish_s_n + 1, s_n, s_n, speed_time_c))
    return True, s_n


def save_yamanouchi_matrix(s_n: int, yd: list, ix: tuple, ym: np.ndarray, speed_time: int, mode=None):
    """
    Yamanouchi对换矩阵的落盘格式为：
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/ij(i,j).pkl 或
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/in(i,n).pkl     ->
    {
    "file_name": "Sn/[ν_i]/ij(i,j)",
    "data": matrix_ij 或 matrix_in,           # list(list(float))
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    ij表示临近对换;
    in表示末尾对换;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1]/ij(2, 3).pkl: {(2,3): [[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]]}
    S3/[2, 1]/in(1, 3).pkl: {(1,3): [[-0.5, -0.8660254037844386], [-0.8660254037844386, 0.5]]}
    """
    if not isinstance(s_n, int) or s_n <= 1:
        err_msg = "s_n={} with type={} must be int and > 1".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(ix, tuple):
        err_msg = "ix={} with type={} must be tuple".format(ix, type(ix))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(ym, np.ndarray) or not ym.size:
        err_msg = "ym={} with type={} must be real np.ndarray with ym.size={}".format(ym, type(ym), ym.size)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(speed_time, int) or speed_time < 0:
        err_msg = "speed_time={} with type={} must be int and >= 0".format(speed_time, type(speed_time))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(mode, str):
        err_msg = "mode={} with type={} must be str".format(mode, type(mode))
        logger.error(err_msg)
        return False, err_msg

    db_info = YamanouchiMatrixInfo(s_n)
    _, file_name = get_yamanouchi_matrix_file_name(s_n, yd, ix, mode=mode)
    table = {"file_name": file_name,
             "data": ym,           # list(list(float))
             "flags": {"speed_time": speed_time}}
    flag, msg = db_info.insert(table)
    if not flag:
        return flag, msg
    flag, msg = db_info.insert_txt(table, point_key="data")
    if not flag:
        return flag, msg

    return True, None


def save_yamanouchi_matrix_finish_s_n(s_n: int, s_n_speed_time: int, is_check_add_one=False):
    """finish_s_n都存txt副本用来展示"""
    if not isinstance(s_n, int) or s_n <= 1:
        err_msg = "s_n={} with type={} must be int and > 1".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n_speed_time, int) or s_n_speed_time < 0:
        err_msg = "s_n_speed_time={} with type={} must be int and >= 0".format(s_n_speed_time, type(s_n_speed_time))
        logger.error(err_msg)
        return False, err_msg

    flag, finish_s_n_before = get_yamanouchi_matrix_finish_s_n()
    if not flag:
        return flag, finish_s_n_before

    if is_check_add_one:
        if s_n - finish_s_n_before != 1:
            err_msg = "is_check_add_one=True require s_n={} - finish_s_n_before={} == 1".format(s_n, finish_s_n_before)
            logger.error(err_msg)
            return False, err_msg

    db_info = YamanouchiMatrixInfo(0)
    _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
    table = {"file_name": finish_file_name,
             "data": [],
             "flags": {"finish_s_n": s_n,
                       "history_times": {
                           "S{}".format(s_n): s_n_speed_time
                       }}}
    if finish_s_n_before == 1:
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


def get_yamanouchi_matrix_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n=1
    """
    _, finish_file_name = get_yamanouchi_matrix_finish_s_n_name()
    flag, data = YamanouchiMatrixInfo(0).query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get yamanouchi_matrix finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, 1
    finish_s_n = data.get("flags", {}).get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= 1:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and > 0, with data={}".format(finish_s_n, data)
        return False, err_msg


def _calc_s_b_to_s_n_part(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, yt_all_s_n):
    """
    方法：
    利用公式
    a，(i-1,i)|Ym> = +- |Ym>， 当i-1和i同行/列
    b，(i-1,i)|Ym> = 1/σ |Ym> + sqrt(σ^2 - 1)/abs(σ) |Yn>， 当不同行列
    其中，σ表示轴距；|Ym>表示标准基；杨盘是Ym
    步骤：
    1，前置必要参数
    2，计算同行列情况
    3，计算异行列情况
    3.1，计算 1/σ |Ym> 部分
    3.2，计算 sqrt(σ^2 - 1)/abs(σ) |Yn> 部分
    4，a式中取相邻对换非对角矩阵元永远大于等于0，为Young-Yamanouchi标准相问
    注意：
    1，这里仅仅计算了(n-1,n)这一个对换矩阵，它既是(ij)又是(in)
    2，轴距按定义，一定是从后面的杨盘（|Ym>），的指标i-1到达i，沿着杨盘向上或向右移动一格+1，沿着杨盘向下或向左移动一格-1，反过来符号就错了！
    """
    # 必须接calc_yamanouchi_matrix_iter，才能不检查入参
    # 1，前置必要参数
    matrix_np = np.zeros([matrix_div, matrix_div])

    # 外层是按照当前young_diagram分支律的分支数bl_num开启循环，其中，range(bl_num)与before_yd一一对应
    # 内层是按照当前分支before_yt_i_batch_len开启循环，循环节是s_b的单个杨盘py序，对Sn的意义是block内部的序号
    # 循环中i标志对于Sn的、整体的序号；j标志对于Sb的、块内的序号
    total_num = 0  # 过程中表示每个分快开始位置
    for bl_index, before_yd_i in zip(range(bl_num), before_yd):
        flag, before_yt_i_batch_len = load_young_table_num(s_b, before_yd_i, is_flag_true_if_not_s_n=False)
        if not flag:  # 这里简化检查手续了
            err_msg = "get before_yt_i_batch db for calc yamanouchi_matrix " \
                      "by s_b={}, before_yt_i_batch_len={} meet error " \
                      "with msg={}".format(s_b, before_yd_i, before_yt_i_batch_len)
            logger.error(err_msg)
            return False, err_msg
        for block_j in range(before_yt_i_batch_len):  # yt yt yt not yd yd yd!
            math_block_j = block_j + 1
            total_i = total_num + block_j
            math_total_i = total_num + math_block_j
            yt_i = yt_all_s_n.get(str(math_total_i))
            row_i_s_b, col_i_s_b = get_s_i_index_in_young_table(s_b, yt_i)
            row_i_s_n, col_i_s_n = get_s_i_index_in_young_table(s_n, yt_i)

            # 2，计算同行列情况
            if row_i_s_b == row_i_s_n:  # 同行最先，且|Yn>为零
                matrix_np[total_i, total_i] = 1
            elif col_i_s_b == col_i_s_n:  # 同列次先，且|Yn>为零
                matrix_np[total_i, total_i] = -1

            else:  # 不同行也不同列，(i-1,i)|Ym> = 1/σ |Ym> + sqrt(σ^2 - 1)/abs(σ) |Yn> 两项都不为零
                # 3，计算异行列情况
                # 3.1，计算 1/σ |Ym> 部分
                axial_distance = row_i_s_b - col_i_s_b + col_i_s_n - row_i_s_n
                matrix_np[total_i, total_i] = 1/axial_distance

                # 3.2，计算 sqrt(σ^2 - 1)/abs(σ) |Yn> 部分  # 当且仅当else分支不为零
                yt_i_after_exchange = copy.deepcopy(yt_i)
                yt_i_after_exchange[row_i_s_n][col_i_s_n] = s_b
                yt_i_after_exchange[row_i_s_b][col_i_s_b] = s_n
                # TODO 在未证明异行列情况下Sb和Sn交换后必然生成合法杨图/Yamanouchi基之前，我们先判断后计算
                # math_another_yt_i_idp_order就是|Yn>的编号
                flag, math_another_yt_i_idp_order = read_young_table_in_decreasing_page_order(
                    yt_i_after_exchange, yt_all_s_n, is_rst_int=True)
                if not flag:
                    err_msg = "read_young_table_in_decreasing_page_order meet error " \
                              "by yt=yt_i_after_exchange={}, yt_all_s_n{} " \
                              "with msg={}".format(yt_i_after_exchange, yt_all_s_n, math_another_yt_i_idp_order)
                    logger.error(err_msg)
                    return False, err_msg
                if math_another_yt_i_idp_order is None:
                    flag, math_another_yt_i_idp_order = quickly_calc_young_table_in_decreasing_page_order(
                        yt_i_after_exchange, yd=before_yd_i, s_n=s_n, is_rst_int=True)
                    if not flag:
                        err_msg = "quickly_calc_young_table_in_decreasing_page_order meet error " \
                                  "by yt=yt_i_after_exchange={}, yd=before_yd_i={}, s_n=s_n={} " \
                                  "with msg={}".format(yt_i_after_exchange, before_yd_i, s_n,
                                                       math_another_yt_i_idp_order)
                        logger.error(err_msg)
                        return False, err_msg
                another_yt_i_idp_order = math_another_yt_i_idp_order - 1
                # sqrt(x^2 - 1)/|x| = sqrt(1 - 1/x^2)
                matrix_np[total_i, another_yt_i_idp_order] = math.sqrt(1 - matrix_np[total_i, total_i] ** 2)

        total_num += before_yt_i_batch_len  # 把结束的block内杨盘总数加进去

    # double check
    if total_num != matrix_div:  # 这个判断是一个分块和等于总数的安全性检查
        err_msg = "total_num={} must eq matrix_div={}, pls check, with s_b={}, s_n={}, yd_s_n={}".format(
            total_num, matrix_div, s_b, s_n, yd_s_n)
        logger.error(err_msg)
        return False, err_msg

    # 注意：这里仅仅计算了(n-1,n)这一个对换矩阵，它既是(ij)又是(in)
    return True, matrix_np


def _calc_remain_matrix_ij_single(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, key_ij):
    """利用分支律，直和出其他(ij)对换矩阵"""
    # 必须接calc_yamanouchi_matrix_iter，才能不检查入参
    # 前置必要参数
    matrix_ij = np.zeros([matrix_div, matrix_div])

    total_num = 0
    for bl_index, before_yd_i in zip(range(bl_num), before_yd):
        flag, matrix_ij_block_i = load_yamanouchi_matrix(s_b, before_yd_i, key_ij,
                                                         mode="ij", is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "load_yamanouchi_matrix meet error, " \
                      "with s_b={}, before_yd_i={}, key_ij={}, mode='ij', " \
                      "msg={}".format(s_b, before_yd_i, key_ij, matrix_ij_block_i)
            logger.error(err_msg)
            return False, err_msg
        flag, matrix_div_ij_block_i = load_young_table_num(s_b, before_yd_i, is_flag_true_if_not_s_n=False)
        if not flag:
            err_msg = "load_young_table_num meet error, with s_b={}, before_yd_i={}, " \
                      "msg={}".format(s_b, before_yd_i, matrix_div_ij_block_i)
            logger.error(err_msg)
            return False, err_msg
        if not matrix_ij_block_i.shape[0] == matrix_ij_block_i.shape[1] == matrix_div_ij_block_i:
            err_msg = "matrix_ij_block_i.shape={} must all eq matrix_div_ij_block_i={} but not".format(
                matrix_ij_block_i.shape, matrix_div_ij_block_i)
            logger.error(err_msg)
            return False, err_msg
        row_tail = total_num + matrix_div_ij_block_i
        matrix_ij[total_num:row_tail, total_num:row_tail] = matrix_ij_block_i
        total_num = row_tail

    # double check
    if total_num != matrix_div:  # 这个判断是一个分块和等于总数的安全性检查
        err_msg = "total_num={} must eq matrix_div={}, pls check, with s_b={}, s_n={}, yd_s_n={}".format(
            total_num, matrix_div, s_b, s_n, yd_s_n)
        logger.error(err_msg)
        return False, err_msg

    return True, matrix_ij


def calc_yamanouchi_matrix_iter(s_n: int, young_diagram: list):
    """
    计算Yamanouchi(ij)(in)对换矩阵
    输入：
    s_n==============>表示当前置换群阶数(int)           # 3
    young_diagram====>表示一个杨图ri(list(int))        # [2, 1]
    输出：(iter)
    y_matrix_ij======>杨图ri的(ij)交换矩阵(matrix)     #{(2,3): [[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]]}
    # 顺序：(Sn-1,Sn) (Sn-2,Sn-1) ... (1,2)
    y_matrix_in======>杨图ri的(in)交换矩阵(matrix)     #{(1,3): [[-0.5, -0.8660254037844386], [-0.8660254037844386, 0.5]]}
    # 顺序：(Sn-1,Sn) (Sn-2,Sn) ... (1,Sn)
    返回：
    flag, {(i,j): y_matrix_ij}/msg, {(i,n): y_matrix_in}/ex_msg
    算法：
    1，先按照定义，赋值(Sn-1，Sn)对换矩阵（适当利用分支律可以优化计算）
    2，利用分支律，直和出其他所有(ij)对换矩阵
    3，使用(ij)计算对应(in)交换矩阵
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        yield False, err_msg, None
        return  # 重点是：下一次迭代时，从上一次迭代遇到的yield后面的代码(下一行)开始执行。
    if not isinstance(young_diagram, list):
        err_msg = "young_diagram={} with type={} must be list".format(young_diagram, type(young_diagram))
        logger.error(err_msg)
        yield False, err_msg, None
        return
    if s_n != sum(young_diagram):
        err_msg = "sum_young_diagram={} must eq s_n={} but not, with young_diagram={}".format(
            sum(young_diagram), s_n, young_diagram)
        logger.error(err_msg)
        yield False, err_msg, None
        return

    # 保护
    if s_n == 1:  # 这个S1在create处也应该有跳过保护
        yield True, False, "s_n=1 do not have Yamanouchi Matrix"
        return
    if s_n == 2:
        if young_diagram == [2]:
            ij_dict = {(1, 2,): np.array([[1]])}
            in_dict = {(1, 2,): np.array([[1]])}
            yield True, ij_dict, in_dict
            return
        elif young_diagram == [1, 1]:
            ij_dict = {(1, 2,): np.array([[-1]])}
            in_dict = {(1, 2,): np.array([[-1]])}
            yield True, ij_dict, in_dict
            return
        else:
            err_msg = "s_n=2 only support yd [2] or [1, 1], but input yd={}".format(s_n, young_diagram)
            logger.error(err_msg)
            yield False, err_msg, None
            return

    # 0，前置数据
    # 分支律数据
    flag, bl_tuple = load_branching_law(s_n, young_diagram, is_flag_true_if_not_s_n=False)
    if not flag:  # 这里简化检查手续了
        err_msg = "get branching_law db for calc yamanouchi_matrix meet error with s_n={}, msg={}".format(s_n, bl_tuple)
        logger.error(err_msg)
        yield False, err_msg, None
        return
    bl_num, rows, cols, before_yd = bl_tuple
    s_b = s_n - 1

    # 杨图杨盘数据
    yd_s_n = copy.deepcopy(young_diagram)
    flag, yt_all_s_n = load_young_table(s_n, yd_s_n, is_flag_true_if_not_s_n=False)
    if not flag:  # 这里简化检查手续了
        err_msg = "get young_table db for calc yamanouchi_matrix meet error with s_n={}, msg={}".format(s_n, yt_all_s_n)
        logger.error(err_msg)
        yield False, err_msg, None
        return

    # 计算矩阵维度矩阵
    flag, matrix_div = load_young_table_num(s_n, yd_s_n, is_flag_true_if_not_s_n=False)
    if not flag:
        err_msg = "load_young_table_num meet error, with s_n={}, yd_s_n={}, msg={}".format(s_n, yd_s_n, matrix_div)
        logger.error(err_msg)
        yield False, err_msg, None
        return

    # 注意，无论(ij)还是(in)，都是二循环（对换），都是群的IR基，指的都是对应 格子 格子 格子 的交换！不是杨图编号的交换！

    # 真正开始计算对换矩阵
    # 1，先按照定义，赋值(Sn - 1，Sn)对换矩阵（适当利用分支律可以优化计算）
    key_ij = key_in = (s_b, s_n,)
    flag, matrix_in = _calc_s_b_to_s_n_part(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, yt_all_s_n)
    if not flag:  # 这里简化检查手续了
        err_msg = "_calc_s_b_to_s_n_part meet error with matrix_div={}, s_b={}, s_n={}, yd_s_n={}, bl_num={}, " \
                  "before_yd={}, yt_all_s_n={}, msg={}".format(matrix_div, s_b, s_n, yd_s_n, bl_num,
                                                               before_yd, yt_all_s_n, matrix_in)
        logger.error(err_msg)
        yield False, err_msg, None
        return

    # 注意：这里仅仅计算了(n-1,n)这一个对换矩阵，它既是(ij)又是(in)
    yield True, {key_ij: matrix_in}, {key_in: matrix_in}  # bingo 1!

    # 注意：不要丢弃这时的matrix_in，它是计算后面其他in的被乘数

    # 建立key头，本身就是倒序哦
    remain_key_ij_list = [(s_b - i - 1, s_n - i - 1,) for i in range(s_n - 2)]
    remain_key_in_list = [(s_b - i - 1, s_n,) for i in range(s_n - 2)]

    for key_ij, key_in in zip(remain_key_ij_list, remain_key_in_list):

        # 2，利用分支律，直和出其他所有(ij)对换矩阵
        flag, matrix_ij = _calc_remain_matrix_ij_single(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, key_ij)
        if not flag:
            err_msg = "_calc_remain_matrix_ij_single meet error, " \
                      "by matrix_div={}, s_b={}, s_n={}, yd_s_n={}, bl_num={}, before_yd={}, key_ij={}" \
                      "with msg={}".format(matrix_div, s_b, s_n, yd_s_n, bl_num, before_yd, key_ij, matrix_ij)
            logger.error(err_msg)
            yield False, err_msg, None
            return

        # 3，使用(ij)计算对应(in)交换矩阵
        # 利用公式 (i, i+v) = (i+1, i+v)(i, i+1)(i+1, i+v)
        # 令n=i+v, j=i+1,有
        # (i, n) = (j, n)(i, j)(j, n)
        matrix_jn = matrix_in  # 式中=右边的in是上一轮的，对于本轮，便是jn
        matrix_in = matrix_jn @ matrix_ij @ matrix_jn  # ij是本轮的
        '''
        The matmul function implements the semantics of the @ operator introduced in Python 3.5 following PEP 465.
        
        matmul differs from dot in two important ways:
            a, Multiplication by scalars is not allowed, use * instead.
            b, Stacks of matrices are broadcast together as if the matrices were elements, 
               respecting the signature (n,k),(k,m)->(n,m)
        
        https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul
        '''

        yield True, {key_ij: matrix_ij}, {key_in: matrix_in}  # bingo 2!

    return


def load_yamanouchi_matrix(s_n: int, yd: list, ix: tuple, mode=None, is_flag_true_if_not_s_n=True):
    """
    取得s_n下young_diagram (ix)的 Yamanouchi矩阵
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    mode_list = ["ij", "in"]  # 这里涉及到调用者必须清楚物理意义，所以不设置"auto"
    if mode not in mode_list:
        err_msg = "mode={} only supported in mode_list={}".format(mode, mode_list)
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(yd, list):
        err_msg = "yd={} with type={} must be list".format(yd, type(yd))
        logger.error(err_msg)
        return False, err_msg
    if not isinstance(ix, tuple):
        err_msg = "ix={} with type={} must be tuple".format(ix, type(ix))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_yamanouchi_matrix_file_name(s_n, yd, ix, mode=mode)
    if not flag:
        err_msg = "cannot get file_name by s_n={}, yd={}, ix={}, mode={} because {}".format(
            s_n, yd, ix, mode, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = YamanouchiMatrixInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query young table with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        yamanouchi_matrix = data.get("data")
        if yamanouchi_matrix.size > 0 and isinstance(yamanouchi_matrix, np.ndarray):
            # 只检查有没有 不对内容做检查了
            return True, yamanouchi_matrix  # bingo！

        else:
            err_msg = "yamanouchi_matrix queried from db, but cannot get yamanouchi_matrix from data" \
                      "with data={}, yamanouchi_matrix={} from db".format(data, yamanouchi_matrix)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_flag_true_if_not_s_n:
            return True, False
        else:
            err_msg = "query not exist yamanouchi_matrix db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg
