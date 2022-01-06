# -*- coding: utf-8 -*-
"""
返回设置好的logger object
TODO 后面需要写多进程并行的话，引入multiprocessing_logging模块
TODO 也可以把log注入StreamHandler，主动触发的时候才写log，手动测试或者当作模块被外部程序引用的时候不写
"""


# import colorlog
import time
import os
import logging
from logging import config  # 这个地方有点没明白，为什么logging不把config暴露出来，难道是不推荐这么用么(python3.6.x)
from functools import lru_cache
from conf.logging_config import logging_config

log_name = logging_config.get("handlers", {}).get("file_handler", {}).get("filename", None)
abs_path = os.path.abspath(__file__)  # /Users/guiyu/PycharmProjects/CGC_of_Sn/src/utils/log.py
log_name = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(abs_path))), log_name)
if log_name.endswith("txt"):
    # hard code to join time mark in log filename
    log_start_time = time.asctime().replace(" ", "_")
    logging_config["handlers"]["file_handler"]["filename"] = \
        log_name.replace("txt", "txt.{}".format(log_start_time))
logging.config.dictConfig(logging_config)


@lru_cache(maxsize=1024)
def get_logger(module_name):
    try:
        logger = logging.getLogger(module_name)
    except Exception as e:
        print(e)
        raise
    return logger
