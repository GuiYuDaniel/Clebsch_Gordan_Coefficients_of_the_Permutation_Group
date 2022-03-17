# -*- coding: utf-8 -*-
"""
this code for creating Branching lengthaw by young_diagrams
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


import copy
from conf.cgc_config import default_s_n
from utils.log import get_logger


logger = get_logger(__name__)


def create_branching_laws(s_n=default_s_n):
    """
    提供给workflow的函数，负责调用计算和存储分支律实体
    返回格式：
    flag, msg
    1，合法：True, s_n
    2，非法：False, msg
    """
    pass


def save_branching_laws():
    pass


def save_branching_laws_finish_s_n():
    pass


def get_branching_laws_finish_s_n():
    pass


def get_branching_laws():
    pass


def calc_single_branching_law(young_diagram: list):
    """
    计算分支律
    输入：
    young_diagram====>表示一个杨图ri(list(int))        # [2, 2]
    输出：
    bl_num===========>表示分支个数(int)                # 1
    rows=============>方格所在行的列表(list(int))       # [1]
    cols=============>方格所在列的列表(list(int))       # [1]
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
        r = copy.deepcopy(young_diagram)
        bl_num = 1
        rows = []
        cols = []
        before_yd = []
        length = len(r)
        if length < 1:
            err_msg = "len={} of r={} must >= 1".format(length, r)
            logger.error(err_msg)
            return False, err_msg
    
        # 计算倒数第一行的分支律
        last_row = length - 1
        if not isinstance(r[last_row], int) or r[last_row] < 1:
            err_msg = "r[{}]={} should be int and >= 1, with r={}".format(last_row, r[last_row], r)
            logger.error(err_msg)
            return False, err_msg
        if r[last_row] == 1:
            rows += [last_row]
            cols += [0]
            before_yd += [r[:-1]]
        else:
            r[last_row] += -1
            rows += [last_row]
            cols += [r[last_row]]
            before_yd += [r[:]]
            r[last_row] += 1

        # 计算倒数第二行到正数第一行的分支律
        if length >= 2:  # 有倒数第二行才有计算
            for i in range(length - 2, -1, -1):  # 从倒数第二行直到第一行 # [l-2, l-3, ..., 0]
                if not isinstance(r[i], int) or r[i] < 1:
                    err_msg = "r[{}]={} should be int and >= 1, with r={}".format(i, r[i], r)
                    logger.error(err_msg)
                    return False, err_msg
                if r[i] > r[i + 1]:  # 本行比下面一行多，才能合法去掉一个格子
                    r[i] += -1  # 真去掉一个，方便后面操作
                    rows += [i]
                    cols += [r[i]]
                    before_yd += [r[:]]
                    bl_num += 1  # 计数器加一
                    r[i] += 1  # 还原成原始的r，再开始下次循环
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
