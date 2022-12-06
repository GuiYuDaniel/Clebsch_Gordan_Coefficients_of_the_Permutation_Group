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
from logging import config
from functools import lru_cache
from conf.logging_config import logging_config
from utils.config import singleton_config

log_name = logging_config.get("handlers", {}).get("file_handler", {}).get("filename", None)  # "logs/log.txt"
log_folder = singleton_config.log_folder
if not os.path.exists(log_folder):
    raise FileExistsError("can not find log_folder={}, pls check".format(log_folder))
log_name = os.path.join(log_folder, log_name)  # p.s. /Users/guiyu/Workspace/CGC/LOG/logs/log.txt
if not os.path.exists(os.path.dirname(log_name)):
    os.makedirs(os.path.dirname(log_name))
if log_name.endswith("txt"):
    # hard code to join time mark in log filename
    log_start_time = time.asctime().replace(" ", "_")
    logging_config["handlers"]["file_handler"]["filename"] = \
        log_name.replace("txt", "txt.{}".format(log_start_time))

if singleton_config.log_level not in ["DEBUG", "INFO"]:
    raise NameError("log_level={} must in ['DEBUG', 'INFO'] but not, pls check config.json"
                    "".format(singleton_config.log_level))
logging_config["loggers"][""]["level"] = singleton_config.log_level

logging.config.dictConfig(logging_config)


@lru_cache(maxsize=1024)
def get_logger(module_name):
    try:
        logger = logging.getLogger(module_name)
    except Exception as e:
        print(e)
        raise
    return logger
