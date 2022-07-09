# -*- coding: utf-8 -*-
"""
python module logging 对应的配置文件，使用py格式而不是json，是为了方便传递
TODO 未来把log的落盘配置写在config里，只留一个需修改文件
"""


logging_config = {
    "version": 1,
    "disable_existing_loggers": True,
    "formatters": {
        "simple": {
            "()": "colorlog.ColoredFormatter",
            "format": "%(log_color)s [%(levelname)s] %(message)s"
        },
        "verbose": {
            "()": "colorlog.ColoredFormatter",
            "format": "%(log_color)s <%(asctime)s> [%(levelname)s] %(message)s [%(filename)s:%(funcName)s:%(lineno)d]",
            "datefmt": "%Y-%m-%d %H:%M:%S"
        }
    },
    "handlers": {
        "console": {
            "level": "INFO",  # 控制台输出的log级别，支持修改为[ERROR, WARNING, INFO, DEBUG]
            "class": "logging.StreamHandler",
            "formatter": "verbose"
        },
        "file_handler": {
            "level": "DEBUG",  # log文件记录的log级别，支持修改为[ERROR, WARNING, INFO, DEBUG]
            "class": "logging.FileHandler",
            "formatter": "verbose",
            "filename": "logs/log.txt",  # log文件在<Clebsch_Gordan_Coefficients_of_the_Permutation_Group>目录下的相对位置，
            # 仅支持改写"logs/log"，后面的".txt"不可修改，
            # 且改写部分不可以 不可以 不可以出现"txt"字样！
            "encoding": "utf8"
        }
    },
    "loggers": {
        "": {
            "handlers": ["console", "file_handler"],
            "level": "INFO"
            # "level": "DEBUG"
        }
    }
}
