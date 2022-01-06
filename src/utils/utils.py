# -*- coding:utf8 -*-
"""
一些其他工具

当需要细化或者函数太多时，应该把其中一些独立出去
"""


import uuid
from enum import Enum, unique
from utils.log import get_logger


logger = get_logger(__name__)


@unique
class PipeTaskStatus(Enum):  # 还不够需要写状态机，先用这个凑活一下
    """
    name is Status
    value is Next
    """
    # TODO 多线程时应该增加STOPPING, WAITING等状态
    # TODO 为了避免Enum键值对别名的规则，引入一个None占位，None无意义，不是一个状态！
    PREPARATION = ["DOING", "FAIL"]  # preparation是准备阶段，计算、保存一切运行的必要信息，未准备结束的ppt是不支持restart的
    DOING = ["SUCCESS", "FAIL", "RESTARTING"]
    SUCCESS = []
    FAIL = ["RESTARTING"]
    RESTARTING = ["DOING", "FAIL", None]  # restarting是重启的准备阶段，准备好后直接进入doing
    # STOPPING = []
    # WAITING = []


def new_id(is_log=True):
    _id = str(uuid.uuid4())
    if is_log:
        logger.debug("new an id={} with uuid4 method".format(_id))
    return _id
