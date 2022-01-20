# -*- coding: utf-8 -*-


import os
import pickle
from functools import lru_cache
from utils.log import get_logger


logger = get_logger(__name__)


class Path(object):

    @classmethod
    def _get_full_path(cls, relative_path="", base_path_type="top"):
        """
        不开放这段代码了，将top，src，test的dirpath登记到Config里，使用时从conf拿更合适
        join the base_path and the relative_path, which base_path was auto fond
        不检查结果，因为外部使用exist调用很方便
        :param relative_path: should be a relative folder or file path(str)
        such as conf, conf/config.json
        :param base_path_type: [top, src, test](str)
        top means <CGC_of_Sn>,
        src means <CGC_of_Sn>/src,
        test means <CGC_of_Sn>/src
        :return: full path(str)
        """
        base_path_type_dict = {"top": "", "src": "src", "test": "test"}
        if not isinstance(relative_path, str):
            err_msg = "relative_path={} should be str".format(relative_path)
            logger.error(err_msg)
            raise Exception(err_msg)
        if base_path_type not in list(base_path_type_dict.keys()):
            err_msg = "base_path_type={} should be in {}".format(base_path_type, list(base_path_type_dict.keys()))
            logger.error(err_msg)
            raise Exception(err_msg)
        base_path = cls._get_top_path()
        full_path = os.path.join(base_path, base_path_type_dict.get(base_path_type), relative_path)
        return full_path

    @staticmethod
    @lru_cache()
    def _get_top_path():
        """
        guess from abspath
        :return: the path of <CGC_of_Sn>(str)
        """
        abs_path = os.path.abspath(__file__)  # <CGC_of_Sn>/src/utils/io.py
        if "src/utils/io.py" in abs_path:
            top_path = os.path.dirname(os.path.dirname(os.path.dirname(abs_path)))
            logger.debug("top path is {} lru cached".format(top_path))
        else:
            logger.warning('Cannot get top_path from abs_path={}'.format(abs_path))
            top_path = None
        return top_path

    def _check_top_path(self, python_path, abs_path):
        """

        :return:
        """
        pass


class Save(object):
    """
    save data
    """

    @staticmethod
    def save_pickle(data, file_path):
        # 范围最宽的path检查应由abc.ABC完成，支持str, bytes or os.PathLike object。我们这里仅允许str就够了
        if not isinstance(file_path, str):
            err_msg = "file_path={} with type={} must be str when save_pickle".format(file_path, type(file_path))
            logger.error(err_msg)
            return False, err_msg
        folder_path = os.path.dirname(file_path)
        if not os.path.exists(folder_path):
            logger.debug("folder_path={} not exist, create it".format(folder_path))
            os.makedirs(folder_path)
        if os.path.exists(file_path):
            err_msg = "file_path={} existed, will not create again, pls check"
            logger.error(err_msg)
            return False, err_msg
        try:
            with open(file_path, "wb") as f:
                pickle.dump(data, f)
                logger.debug("save data={} to db={}".format(data, file_path))
        except Exception as e:
            logger.error("save data to file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e
        return True, True


class Load(object):
    """
    load data
    """

    @staticmethod
    def load_pickle(file_path):
        # 范围最宽的path检查应由abc.ABC完成，支持str, bytes or os.PathLike object。我们这里仅允许str就够了
        if not isinstance(file_path, str):
            err_msg = "file_path={} with type={} must be str when load_pickle".format(file_path, type(file_path))
            logger.error(err_msg)
            return False, err_msg
        if not os.path.exists(file_path):
            err_msg = "file_path={} not existed, will not load, pls check"
            logger.error(err_msg)
            return False, err_msg
        if not os.path.isfile(file_path) or file_path[-4:] != ".pkl":  # 先用名字查类型
            err_msg = "file_path={} must a pkl file, cannot load, pls check".format(file_path)
            logger.error(err_msg)
            return False, err_msg
        try:
            with open(file_path, "rb") as f:
                data = pickle.load(f)
                logger.debug("load data={} from db={}".format(data, file_path))
        except Exception as e:
            logger.error("load data from file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e
        return True, data


class Delete(object):
    """
    delete data
    """

    @staticmethod
    def delete_pickle(file_path):
        # 范围最宽的path检查应由abc.ABC完成，支持str, bytes or os.PathLike object。我们这里仅允许str就够了
        if not isinstance(file_path, str):
            err_msg = "file_path={} with type={} must be str when delete_pickle".format(file_path, type(file_path))
            logger.error(err_msg)
            return False, err_msg
        if not os.path.exists(file_path):
            err_msg = "file_path={} not existed, cannot delete but return True"
            logger.error(err_msg)
            return False, err_msg
        try:
            os.remove(file_path)
            logger.debug("delete data from db={}".format(file_path))
        except Exception as e:
            logger.error("delete data from file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e
        return True, True
