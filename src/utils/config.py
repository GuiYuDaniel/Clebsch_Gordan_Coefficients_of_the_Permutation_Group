# -*- coding: utf-8 -*-


import os
from utils.json_op import Json
# from utils.log import get_logger
#
#
# logger = get_logger(__name__)


class Config(object):
    """
    An useful config class, with cmd > user defined > default
    """

    def __init__(self):
        # general config
        self.name = None
        self.version = None
        self.github = None
        self.paper = None
        self.author = None
        self.license = None

        # system config  # user defined
        self.log_folder = None
        self.log_level = None
        self.result_folder = None
        self.sys_db_name = None

        # core config  # user defined
        self.calc_target = None

        # project config
        self.top_path = self._get_top_path()
        self.src_path = os.path.join(self.top_path, "src")
        self.test_path = os.path.join(self.top_path, "test")

        default_config = self._load_default_config()
        self._add_default_config(default_config)
        self.add_cmd_config()
        self._last_check()

    @staticmethod
    def _get_top_path():
        """guess from abspath"""
        abs_path = os.path.abspath(__file__)  # <top_path>/src/utils/config.py
        file_path = os.path.join("src", "utils", "config.py")
        if file_path in abs_path:
            top_path = os.path.dirname(os.path.dirname(os.path.dirname(abs_path)))
            return top_path  # <top_path>
        else:
            raise Exception("Cannot get top_path from abs_path={}".format(abs_path))

    def _load_default_config(self):
        """Not all config params will be added"""
        default_config_path = os.path.join(self.src_path, "conf", "config.json")
        default_config = Json.file_to_json_without_comments(default_config_path)
        return default_config

    def _add_default_config(self, default_config):
        self.name = default_config.get("name")
        self.version = default_config.get("version")
        self.github = default_config.get("github")
        self.paper = default_config.get("paper")
        self.author = default_config.get("author")
        self.license = default_config.get("license")

        self.log_folder = default_config.get("DEFAULT_LOG_FOLDER")
        self.log_level = default_config.get("DEFAULT_LOG_LEVEL")
        self.result_folder = default_config.get("DEFAULT_RESULT_FOLDER")
        self.sys_db_name = default_config.get("DEFAULT_DB_NAME")
        self.calc_target = default_config.get("DEFAULT_CALC_TARGET")

    def add_cmd_config(self):
        pass

    def _last_check(self):
        none_dict = {k: self.__dict__[k] for k in self.__dict__ if self.__dict__[k] is None}
        if none_dict:
            err_msg = "config must have no None but {}".format(none_dict)
            # logger.error(err_msg)
            raise Exception(err_msg)

        # logger.info("Checked config is {}".format(self.__dict__))


singleton_config = Config()
