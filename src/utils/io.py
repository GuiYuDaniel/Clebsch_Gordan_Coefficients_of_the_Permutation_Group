# -*- coding: utf-8 -*-


import os
import pickle
from functools import lru_cache, wraps, partial
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
        test means <CGC_of_Sn>/test
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
        # 除去""的目的是避免join("/xxx", "")返回"/xxx/"这样的结果
        real_partial_path_list = [i for i in [base_path, base_path_type_dict.get(base_path_type), relative_path] if i]
        full_path = os.path.join(*tuple(real_partial_path_list))
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


def _make_file_path_with_lawful_file_type(file_path, file_type):
    # 检查file_type
    if not isinstance(file_type, str) or not file_type.startswith("."):
        err_msg = "file_type={} must a str and startswith '.'".format(file_type)
        logger.debug(err_msg)
        return False, err_msg
    # 范围最宽的path检查应由abc.ABC完成，支持str, bytes or os.PathLike object 我们这里仅允许str就够了
    if not isinstance(file_path, str):
        err_msg = "file_path={} with type={} must be str when save_pickle".format(file_path, type(file_path))
        logger.debug(err_msg)
        return False, err_msg
    # 添加一个abs检查，毕竟return的合格路径马上就要开始使用了，它仅仅检查绝对性，所以和扩展名有无无关
    if not os.path.isabs(file_path):
        err_msg = "file_path={} must a abspath but not".format(file_path)
        logger.debug(err_msg)
        return False, err_msg

    # 这段负责确保file_path是一个以.xxx结尾的合法文件路径
    file_name, ext_name = os.path.splitext(file_path)
    if ext_name == "":  # file_path没有带扩展名，加上file_type返回
        file_path = file_path + file_type
        return True, file_path
    if ext_name == file_type:  # 带了正确的扩展名，原样返回
        return True, file_path
    else:  # 带错扩展名了，false
        err_msg = "file_path={} not match point file_type={}, pls check".format(file_path, file_type)
        logger.debug(err_msg)
        return False, err_msg


def _revise_path_before_io(file_type=None, io_type=None):
    """
    装饰器作用：
    1，检查参数合法性
    2，拼装完整路径
    3，创建外部文件夹
    """
    # 检查装饰器参数
    if file_type is None:
        err_msg = "need real file_type but {}".format(file_type)
        logger.error(err_msg)
        return False, err_msg
    lawful_io_type = ["save", "load", "delete"]
    if io_type not in lawful_io_type:
        err_msg = "io_type={} must in lawful_io_type={} but not".format(io_type, lawful_io_type)
        logger.error(err_msg)
        return False, err_msg

    # 不带参数的装饰器写法
    def revise_path_before_io(func):
        if not callable(func):
            er_msg = "function name={} is not callable, pls check".format(func.__name__)
            logger.error(er_msg)
            return False, er_msg
        wrapper = None, None

        # save的前置逻辑
        if io_type == "save":
            @wraps(func)
            def wrapper(data, file_path):
                flag, file_path = _make_file_path_with_lawful_file_type(file_path, file_type)
                if not flag:
                    e_msg = file_path
                    logger.error(e_msg)
                    return False, e_msg
                folder_path = os.path.dirname(file_path)
                if not os.path.exists(folder_path):
                    logger.debug("folder_path={} not exist, create it".format(folder_path))
                    os.makedirs(folder_path)
                if os.path.exists(file_path):
                    e_msg = "file_path={} existed, will not create again, pls check".format(file_path)
                    logger.error(e_msg)
                    return False, e_msg
                result = func(data, file_path)
                return result

        # load的前置逻辑
        elif io_type == "load":
            @wraps(func)
            def wrapper(file_path):
                flag, file_path = _make_file_path_with_lawful_file_type(file_path, file_type)
                if not flag:
                    e_msg = file_path
                    logger.error(e_msg)
                    return False, e_msg
                if not os.path.exists(file_path):
                    e_msg = "file_path={} not existed, will not load, pls check".format(file_path)
                    logger.error(e_msg)
                    return False, e_msg
                if not os.path.isfile(file_path):  # 先用名字查类型
                    e_msg = "file_path={} must a file, cannot load, pls check".format(file_path)
                    logger.error(e_msg)
                    return False, e_msg
                result = func(file_path)
                return result

        # delete的前置逻辑
        elif io_type == "delete":
            @wraps(func)
            def wrapper(file_path):
                flag, file_path = _make_file_path_with_lawful_file_type(file_path, file_type)
                if not flag:
                    e_msg = file_path
                    logger.error(e_msg)
                    return False, e_msg
                if not os.path.exists(file_path):
                    e_msg = "file_path={} not existed, cannot delete but return True, msg".format(file_path)
                    logger.warning(e_msg)
                    return True, e_msg
                result = func(file_path)
                return result

        # 由于前面检查了io_type，不会真正走到这个分支，形式逻辑上有它更加易读
        else:
            pass

        return wrapper

    return revise_path_before_io


_revise_path_before_save = partial(_revise_path_before_io, io_type="save")
_revise_path_before_load = partial(_revise_path_before_io, io_type="load")
_revise_path_before_delete = partial(_revise_path_before_io, io_type="delete")


# def _revise_path_before_save(*args, **kwargs):
#     _func = partial(_revise_path_before_io, io_type="save")
#     return _func(*args, **kwargs)


class Save(object):
    """
    save data
    """

    @staticmethod
    @_revise_path_before_save(file_type=".pkl")
    def save_pickle(data, file_path):
        try:
            with open(file_path, "wb") as f:
                pickle.dump(data, f)
                logger.debug("save data={} to db={}".format(data, file_path))
        except Exception as e:
            logger.error("save data to file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e

        return True, True

    @staticmethod
    @_revise_path_before_save(file_type=".txt")
    def save_txt(data, file_path):
        if not isinstance(data, str):
            data = str(data)
        try:
            with open(file_path, "w") as f:
                f.write(data)
                logger.debug("save data={} in={}".format(data, file_path))
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
    @_revise_path_before_load(file_type=".pkl")
    def load_pickle(file_path):
        try:
            with open(file_path, "rb") as f:
                data = pickle.load(f)
                logger.debug("load data={} from db={}".format(data, file_path))
        except Exception as e:
            logger.error("load data from file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e

        return True, data

    # TODO 考虑要不要对data加str检查和强制类型转化
    @staticmethod
    @_revise_path_before_load(file_type=".txt")
    def load_txt(file_path):
        try:
            with open(file_path, "r") as f:
                data = f.read()
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
    @_revise_path_before_delete(file_type=".pkl")
    def delete_pickle(file_path):
        try:
            os.remove(file_path)
            logger.debug("delete data from db={}".format(file_path))
        except Exception as e:
            logger.error("delete data from file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e
        return True, True

    @staticmethod
    @_revise_path_before_delete(file_type=".txt")
    def delete_txt(file_path):
        try:
            os.remove(file_path)
            logger.debug("delete data from {}".format(file_path))
        except Exception as e:
            logger.error("delete data from file_path={} meet error".format(file_path))
            logger.error(e)
            return False, e
        return True, True
