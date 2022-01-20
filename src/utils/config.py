# -*- coding: utf-8 -*-


from utils.io import Path
from utils.json_op import Json
from utils.log import get_logger


logger = get_logger(__name__)


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
        # TODO 把log实装到conf里
        # self.log_folder = None
        # self.log_level = None
        self.result_folder = None

        # core config  # user defined
        self.calc_target = None

        # project config
        self.top_path = Path._get_top_path()
        self.src_path = Path._get_full_path(base_path_type="src")
        self.test_path = Path._get_full_path(base_path_type="test")

        default_config = self._load_default_config()
        self._add_default_config(default_config)
        self.add_cmd_config()
        self._last_check()

    @staticmethod
    def _load_default_config():
        """Not all config params will be added"""
        default_config_path = Path._get_full_path(relative_path="conf/config.json", base_path_type="src")
        default_config = Json.file_to_json_without_comments(default_config_path)
        return default_config

    def _add_default_config(self, default_config):
        self.name = default_config.get("name")
        self.version = default_config.get("version")
        self.github = default_config.get("github")
        self.paper = default_config.get("paper")
        self.author = default_config.get("author")
        self.license = default_config.get("license")

        self.result_folder = default_config.get("DEFAULT_RESULT_FOLDER")
        self.calc_target = default_config.get("DEFAULT_CALC_TARGET")

    def add_cmd_config(self):
        pass

    def _last_check(self):
        none_dict = {k: self.__dict__[k] for k in self.__dict__ if self.__dict__[k] is None}
        if none_dict:
            err_msg = "config must have no None but {}".format(none_dict)
            logger.error(err_msg)
            raise Exception(err_msg)

        logger.info("Checked config is {}".format(self.__dict__))


singleton_config = Config()
