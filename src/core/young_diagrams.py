# -*- coding: utf-8 -*-
"""
this code for creating Young Diagrams by Sn
"""


from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name, get_young_diagrams_finish_s_n_name
from conf.cgc_config import default_s_n
from utils.log import get_logger


logger = get_logger(__name__)


def calc_young_diagrams(s_n=default_s_n):
    """
    提供给workflow的函数，负责计算和存储
    返回格式：
    flag, msg
    1，合法：True, None
    2，非法：False, msg
    """
    if not isinstance(s_n, int) or s_n <= 0:
        err_msg = "s_n={} with type={} must be int and > 0".format(s_n, type(s_n))
        logger.error(err_msg)
        return False, err_msg

    # 先查询数据库中完成到S几：如果s_n超过它，给出循环头；如果s_n被计算过了，则给出完成标记（注意不是返回结果）
    flag, finish_s_n = get_young_diagrams_finish_s_n()
    if not flag:
        err_msg = "get young_diagrams finish s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了
        msg = "s_n={} young_diagrams had been calculated, return True, None".format(s_n)
        logger.info(msg)
        return True, None
    else:
        msg = "finish_s_n={}, will calc s_n from={} to {}".format(finish_s_n, finish_s_n + 1, s_n)
        logger.info(msg)

    # 按照从小到大的顺序，逐个计算s_i的杨图并储存
    for s_i in range(finish_s_n + 1, s_n + 1):  # 循环体为[finish_s_n+1, finish_s_n+2, ..., s_n]
        flag, s_i_young_diagrams = calc_single_young_diagrams(s_i)  # 因为逐次向上计算，所以recursion_deep可以限制为1
        if not flag:
            err_msg = "calc s_i_young_diagrams meet error with s_i={}, msg={}".format(s_i, s_i_young_diagrams)
            logger.error(err_msg)
            return False, err_msg
        if not s_i_young_diagrams:  # 因为逐次向上计算，所以recursion_deep可以限制为1
            err_msg = "calc s_i_young_diagrams should not kill by recursion with s_i={}, pls check".format(s_i)
            logger.error(err_msg)
            return False, err_msg
        # 既然没问题了，那就存（别忘了也要更新Finish_Sn）

    # 返回消息
    pass


def save_young_diagrams(s_n: int, young_diagrams: list, speed_time=None):
    """
    杨图的落盘格式为：
    <top_path>/cgc_results/young_diagrams_info/Sn.pkl ->
    {
    "data": [[gamma_i], ...]
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [[gamma_i], ...]表示配分，是二维列表;
    speed_time表示计算用时（秒）

    例如：
    S3.pkl: [[3], [2, 1], [1, 1, 1]]
    S4.pkl: [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    """
    pass


def save_young_diagrams_finish_s_n(s_n: int):


def calc_single_young_diagrams(s_n: int, recursion_deep: int=1):
    """
    杨图的计算程序：
    根据给定的Sn，计算全部杨图的配分，并返回
    （法五：首行数字循环，后面不递归了，而是读取前面的结果）

    recursion_deep表示当前计算阶数领先于业已计算阶数多少不被cut，
    例如0就只能取结果不会真算；1表示可以比算好的大1

    结果举例；
    S3: [[3], [2, 1], [1, 1, 1]]
    S6: [[6], [5, 1], [4, 2], [4, 1, 1], [3, 3], [3, 2, 1], [3, 1, 1, 1],
         [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]
    注意：排序按照金牌模式而不是行数，例如[3, 3]两行排在了[4, 1, 1]后面
    """
    if not isinstance(s_n, int) or s_n <= 0 or not isinstance(recursion_deep, int) or recursion_deep < 0:
        err_msg = "s_n={} with type={} must be int and > 0 recursion_deep={} with type={} must be int and >= 0".format(
            s_n, type(s_n), recursion_deep, type(recursion_deep))
        logger.error(err_msg)
        return False, err_msg

    # 先查询数据库中完成到S几，如果没计算到，根据recursion_deep判断是否计算；如果计算到了，就必须拿到结果
    flag, finish_s_n = get_young_diagrams_finish_s_n()
    if not flag:
        err_msg = "get young_diagrams finish s_n meet error with msg={}".format(finish_s_n)
        logger.error(err_msg)
        return False, err_msg
    if s_n <= finish_s_n:
        # 说明以前算过了，直接返回结果，拿不到要抛错！
        flag, young_diagrams = load_young_diagrams(s_n, is_return_false_if_not_s_n=False)
        if not flag:
            err_msg = "get young_diagrams db meet error with s_n={}".format(s_n)
            logger.error(err_msg)
            return False, err_msg
        return True, young_diagrams
    elif s_n - finish_s_n <= recursion_deep:
        # 落在recursion_deep范围内，开始下面的递归计算
        pass
    else:
        # 落在recursion_deep范围外，礼貌推出
        logger.debug("s_n={}, finish_s_n={}, deep need {} out of recursion_deep={}, will return True, False".format(
            s_n, finish_s_n, s_n - finish_s_n, recursion_deep))
        return True, False

    # 保护
    if s_n == 1:
        young_diagrams = [[1]]
        return True, young_diagrams

    # 真正开始计算young_diagrams
    young_diagrams = []
    for first_num in range(s_n, 0, -1):  # first_num 即杨图的首行格子数，循环体为[s_n, s_n-1, ..., 1]
        remain_num = s_n - first_num  # 除首行外，剩余的格子数。用它作为Sx来查询前面以后结果，根据first_num直接接上
        young_diagrams_batch = []  # 因为杨图是一波一波算的，每一波叫一个batch

        if first_num == s_n:  # 首行即是全部格子的，那个batch杨图就是它自己。remain_num也必然是0
            young_diagram_single = [first_num]
            young_diagrams_batch.append(young_diagram_single)
            young_diagrams += young_diagrams_batch
            continue

        # 取得除第一层外，剩余阶数的子杨图：sub_young_diagrams
        flag, sub_young_diagrams = calc_single_young_diagrams(remain_num, recursion_deep=recursion_deep)
        if not flag:
            err_msg = "get sub_young_diagrams db meet error with remain_num={}".format(remain_num)
            logger.error(err_msg)
            return False, err_msg
        if not sub_young_diagrams:  #
            err_msg = "calc sub_young_diagrams should recursion with remain_num={}, " \
                      "and recursion_deep for sub must < for itself," \
                      "sub_young_diagrams={} not expect False, pls check".format(remain_num, sub_young_diagrams)
            logger.error(err_msg)
            return False, err_msg

        # 切除子杨图中，首行格子数大于first_num的那些，因为形如[2, 5, 1]这种图不符合杨图规则
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
            err_msg = "cut sub_young_diagrams meet error, " \
                      "with sub_young_diagrams={}, first_num={}".format(sub_young_diagrams, first_num)
            logger.error(err_msg)
            logger.error(e)
            return False, e

        young_diagrams += young_diagrams_batch

    return True, young_diagrams


def get_young_diagrams_finish_s_n():
    """
    flag表示是否有报错，
    finish_s_n表示当前计算完成的Sn，如果没有，则finish_s_n=0
    """
    finish_file_name = get_young_diagrams_finish_s_n_name()
    flag, data = YoungDiagramInfo().query_by_file_name(finish_file_name)
    if not flag:
        err_msg = "get young_diagrams finish_s_n meet error with finish_file_name={}".format(finish_file_name)
        logger.error(err_msg)
        return False, err_msg
    if data is False:
        # logger.debug("find no finish_s_n, return 0")
        return True, 0
    finish_s_n = data.get("flags").get("finish_s_n")
    if finish_s_n and isinstance(finish_s_n, int) and finish_s_n >= 1:
        return True, finish_s_n
    else:
        err_msg = "finish_s_n={} must int and > 0, with data={}".format(finish_s_n, data)
        return False, err_msg


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
        young_diagrams = data.get("data")
        if young_diagrams:
            # 只检查有没有 不对内容做检查了
            # TODO 谁用谁检查，目前就calc_single_young_diagrams用了，它也检查了
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
            err_msg = "query not exist young_diagrams db with s_n={}, file_name={}, err_msg={}".format(
                s_n, file_name, data)
            return False, err_msg


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
