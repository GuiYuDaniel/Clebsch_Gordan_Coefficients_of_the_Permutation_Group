# -*- coding:utf8 -*-
"""
这不是一个严格意义上的数据库，只是使用存取pickle，统一接口形式，满足最简单的增删改查需要
在project的db基础上，多了txt的保存需求，需要改写一些函数
"""


import os
from conf.cgc_config import (
    top_path, cgc_rst_folder
)
from db.local_db import LocalDb
from utils.log import get_logger

logger = get_logger(__name__)


class CGCDb(LocalDb):

    def __init__(self):
        super(CGCDb, self).__init__()
        self.db_folder = os.path.join(top_path, cgc_rst_folder, "{}", "{}")
        self.design_table_type.update({  # db会自动添加create_time:str和last_write_time:str两项
            "pipetask_id": str,

            "data": None,
            "flags": None
        })

    # TODO 重写一些函数，要考虑深度不定，类型可判断两件事
