# -*- coding: utf-8 -*-
"""
this code for creating Young Diagrams by Sn
"""


from core.cgc_utils.cgc_db_typing import CalculatedTableInfo
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name
from conf.cgc_config import default_s_n
from utils.log import get_logger


logger = get_logger(__name__)


def calc_young_diagrams(s_n=default_s_n):
    pass


def calc_single_young_diagrams(s_n: int, recursion_deep: int=1):
    """
    杨图的计算程序：
    根据给定的Sn，计算全部杨图的配分，并返回
    （法五：首行数字循环，后面不递归了，而是读取前面的结果）

    例如S3: [[3], [2, 1], [1, 1, 1]]
    例如S6: [[6], [5, 1], [4, 2], [4, 1, 1], [3, 3], [3, 2, 1], [3, 1, 1, 1],
            [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]
    注意：排序按照金牌模式而不是行数，例如[3, 3]两行排在了[4, 1, 1]后面
    """
    if not isinstance(s_n, int) or s_n <= 0 or not isinstance(recursion_deep, int):
        err_msg = "s_n={} with type={} must be int and > 0 recursion_deep={} with type={} must be int".format(
            s_n, type(s_n), recursion_deep, type(recursion_deep))
        logger.error(err_msg)
        return False, err_msg

    # 先查询数据库中完成到S几，如果没计算到，根据recursion_deep判断是否计算；如果计算到了，就必须拿到结果

    flag, young_diagrams = load_young_diagrams(s_n, is_return_false_if_not_s_n=True)  # 先查询db，看看是否有结果
    if not flag:
        err_msg = "get young_diagrams db meet error with s_n={}".format(s_n)
        logger.error(err_msg)
        return False, err_msg
    if young_diagrams:  # 说明已经计算过 直接返回结果
        return True, young_diagrams
    if not is_recursion:  # 不递归
        return True, False

    if s_n == 1:
        young_diagrams = [[1]]
        return True, young_diagrams

    # 真正开始计算young_diagrams
    young_diagrams = []
    for first_num in range(s_n, 0, -1):  # first_num 即杨图的首行格子数，循环体为[s_n, s_n-1, ..., 1]
        remain_num = s_n - first_num  # 除首行外，剩余的格子数。用它作为Sx来查询前面以后结果，根据first_num直接接上
        young_diagrams_batch = []  # 因为杨图是一波一波算的，每一波叫一个batch

        if first_num == s_n:  # 首行即是全部格子的，那个batch杨图就是它自己。remain_num也必然是0
            single_young_diagrams = [first_num]
            young_diagrams_batch.append(single_young_diagrams)
            young_diagrams += young_diagrams_batch
            continue

        # TODO 是否需要限制递归深度
        flag, sub_young_diagrams = calc_single_young_diagrams(remain_num, is_recursion=True)
        if not flag:
            err_msg = "get sub_young_diagrams db meet error with remain_num={}".format(remain_num)
            logger.error(err_msg)
            return False, err_msg
        if not sub_young_diagrams:  #
            err_msg = "calc sub_young_diagrams should recursion with remain_num={}, " \
                      "sub_young_diagrams={} not expect False, pls check".format(remain_num, sub_young_diagrams)
            logger.error(err_msg)
            return False, err_msg

        # 切除子杨盘中，首行格子数大于first_num的那些，因为形如[2, 5, 1]这种图不符合杨图规则
        try:  # 由于不能保证取得的子杨图肯定是二维列表，所以使用try
            cut_sub_young_diagrams = [i for i in sub_young_diagrams if i[0] <= first_num]
            young_diagrams_batch = [[first_num] + i for i in cut_sub_young_diagrams]
            # TODO log写在调用函数的地方吧，不然因为递归，会写很多
            # logger.debug("calc young_diagrams_batch with \n"
            #              "s_n={}, first_num={}, remain_num={}, \n"
            #              "so, sub_young_diagrams={}, cut_sub_young_diagrams={}, \n "
            #              "young_diagrams_batch={}, will be add to young_diagrams={}".format(
            #     s_n, first_num, remain_num,
            #     sub_young_diagrams, cut_sub_young_diagrams,
            #     young_diagrams_batch, young_diagrams))
        except Exception as e:
            logger.error("cut sub_young_diagrams={} meet error".format(sub_young_diagrams))
            logger.error(e)
            return False, e

        young_diagrams += young_diagrams_batch

    return True, young_diagrams


def save_young_diagrams():
    """
    杨图的落盘格式为：
    young_diagrams/Sn.type，内容为[[gamma_i], ...]

    例如：
    young_diagrams/S3.pkl: [[3], [2, 1], [1, 1, 1]]
    young_diagrams/S4.txt: [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    """
    pass


def load_young_diagrams(s_n: int, is_return_false_if_not_s_n=True):
    """
    取得s_n的杨图(二维列表)
    如果没有，根据is_return_true_if_not_s_n决定返回True or False
    """
    if not isinstance(s_n, int) or not isinstance(is_return_false_if_not_s_n, bool):
        err_msg = "s_n={} with type={} must be int, is_return_false_if_not_s_n={} with type={} must be bool".format(
            s_n, type(s_n), is_return_false_if_not_s_n, type(is_return_false_if_not_s_n))
        logger.error(err_msg)
        return False, err_msg

    flag, file_name = get_young_diagrams_file_name(s_n)
    if not flag:
        err_msg = "cannot get file_name by s_n={} because {}".format(s_n, file_name)
        logger.error(err_msg)
        return False, err_msg
    flag, data = YoungDiagramInfo(s_n).query_by_file_name(file_name)
    if not flag:
        err_msg = "cannot query young diagrams with s_n={}, file_name={} because {}".format(s_n, file_name, data)
        logger.error(err_msg)
        return False, err_msg

    if data:
        young_diagrams = data.get("data").get("young_diagrams")
        if young_diagrams:
            # 只检查有没有 不对内容做检查了
            return True, young_diagrams  # bingo！

        else:
            err_msg = "young_diagrams queried form db, but cannot get young_diagrams from data" \
                      "with data={}, young_diagrams={} from db".format(data, young_diagrams)
            logger.error(err_msg)
            return False, err_msg
    else:
        if is_return_false_if_not_s_n:
            return True, False
        else:
            return False, "query not exist young_diagrams db with s_n={}, file_name={}".format(s_n, file_name)


"""
备选算法：
法一：
按照行数循环（太多if）
法二：
按照首行数字大小递归（还行，就是循环深度随着n线性增大，计算量平方增大）（好处是无需再做排序）
法三：
排列组合（利用python内置函数itertools，可能 可能 可能计算得快）
法四：
倒序生长（不如直接用分支律呢）
法五：
首行数字循环，后面不递归了，而是读取前面的结果（机智啊！io开销并不大，因为只需读取Sn个，
并且，第二行的，前面必然有结果了，不需要递归了，循环深度===2）
法六：[当前！]
修改版法五，如果没找到之前的，递归寻找，增强了单独计算某个远端Sn的能力
"""

