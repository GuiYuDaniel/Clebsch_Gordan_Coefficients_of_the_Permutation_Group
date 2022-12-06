# -*- coding:utf8 -*-
"""
保护生产db，这里面错了，可以用raise打断进程
"""


import os
from utils.config import singleton_config
from utils.log import get_logger


logger = get_logger(__name__)


class DBProtector(object):
    """
    供需要保护生产db的时候使用
    对于pkl实现的local db，原理是：

    mv mode:
    setup
    1，mv当前 db --> db.protected
    2，新建test db
    teardown
    3，删除test db
    4，mv回来 db.protected --> db

    copy mode:
    setup
    cp db --> db.protected, the new test db is db
    teardown
    1, delete test db
    2, mv db.protected --> db
    """

    def __init__(self, relative_db_path: str, is_copy=False,
                 extension_name: str=".protected", full_path_ahead: str=None):
        if not isinstance(relative_db_path, str):
            err_msg = "relative_src={} must be str".format(relative_db_path)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._relative_db_path = relative_db_path

        if full_path_ahead is None:
            self._full_path_ahead = singleton_config.result_folder
        else:
            if not isinstance(full_path_ahead, str):
                err_msg = "full_path_ahead={} must be str".format(full_path_ahead)
                logger.error(err_msg)
                raise Exception(err_msg)
            self._full_path_ahead = full_path_ahead

        if not isinstance(is_copy, bool):
            err_msg = "is_copy={} must be bool".format(is_copy)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._is_copy = is_copy

        if not isinstance(extension_name, str) or not extension_name.startswith("."):
            err_msg = "extension_name={} must be str and startswith '.'".format(extension_name)
            logger.error(err_msg)
            raise Exception(err_msg)

        self._full_src = os.path.join(self._full_path_ahead, self._relative_db_path)
        self._full_dst = self._full_src + extension_name

        # 添加一个状态位，标志是否已经对生产环境做了保护
        self._protected_status = False
        self._db_exists = None

        self._check_init()

    def _check_init(self):
        """
        检查目录安全
        检查被保护的目录是否存在
        检查临时存放目录不存在
        """
        # if singleton_config.top_path not in self._full_src or singleton_config.top_path == self._full_src:
        #     err_msg = "the db path={} must under top_path={}".format(self._full_src, singleton_config.top_path)
        #     logger.error(err_msg)
        #     raise Exception(err_msg)
        ban_str_list = ["..", "*"]  # 这些字符不可以出现在目录中
        for ban_path in ban_str_list:
            if ban_path in self._full_src or ban_path in self._full_dst:
                err_msg = "the db path={} and db.protected path={} must not have {}".format(
                    self._full_src, self._full_dst, ban_path)
                logger.error(err_msg)
                raise Exception(err_msg)
        if self._full_src.startswith(singleton_config.top_path):
            err_msg = "the db path={} should not startswith {}".format(self._full_src, singleton_config.top_path)
            logger.error(err_msg)
            raise Exception(err_msg)

        # 现在生产目录存在/不存在都支持了
        # 如果存在 转移走保护起来 结束后再转移回去
        # 如果不存在 新建一个该目录 结束后删除
        if os.path.exists(self._full_src):
            self._db_exists = True
        else:
            self._db_exists = False

        if os.path.exists(self._full_dst):
            err_msg = "the protector db path={} must be not exists".format(self._full_dst)
            logger.error(err_msg)
            raise Exception(err_msg)

    def protector_setup(self):
        if self._protected_status is True:
            err_msg = "db={} had been protected as {}, cannot protect again".format(self._full_src, self._full_dst)
            logger.error(err_msg)
            raise Exception(err_msg)

        if self._db_exists:
            if self._is_copy:  # copy mode
                # 1, copy db(src) to db.protected(dst), db become to the new test_db
                os.system("cp -r {} {}".format(self._full_src, self._full_dst))
                msg = "db={} has been protected to {} with copy mode".format(self._full_src, self._full_dst)
                logger.info(msg)
            else:  # mv mode
                # 1, mv db(src) to db.protected(dst)
                os.system("mv {} {}".format(self._full_src, self._full_dst))
                # os.rename(self._full_src, self._full_dst)
                # 2, makedir a new blank test_db(src)
                os.makedirs(self._full_src)
                msg = "db={} has been protected to {} with mv mode".format(self._full_src, self._full_dst)
                logger.info(msg)
        else:
            os.makedirs(self._full_src)
            msg = "db={} is not exists, create it, and will delete in teardown for protecting".format(self._full_src)
            logger.info(msg)

        self._protected_status = True

    def protector_teardown(self):
        if self._protected_status is False:
            err_msg = "protected db={} had been teardown to {}, cannot teardown again".format(
                self._full_dst, self._full_src)
            logger.error(err_msg)
            raise Exception(err_msg)

        if self._db_exists:
            # 1, rm test_db(src)
            os.system("rm -rf {}".format(self._full_src))
            # 2, mv db.protected(dst) to db(src)
            os.system("mv {} {}".format(self._full_dst, self._full_src))
            # os.rename(self._full_dst, self._full_src)
            msg = "protector db={} has been mv back to {}".format(self._full_dst, self._full_src)
            logger.info(msg)
        else:
            os.system("rm -rf {}".format(self._full_src))
            msg = "protector db={} has been delete".format(self._full_src)
            logger.info(msg)

        self._protected_status = False
