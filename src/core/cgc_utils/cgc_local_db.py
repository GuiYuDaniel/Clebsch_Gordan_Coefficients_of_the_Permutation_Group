# -*- coding:utf8 -*-
"""
这不是一个严格意义上的数据库，只是使用存取pickle，统一接口形式，满足最简单的增删改查需要
在继承project的local_db基础上，添加一些txt函数用于展示结果
模仿src/db的形式写的
"""


import copy
import os
import sys
from conf.cgc_config import (
    top_path, cgc_rst_folder
)
from db.local_db import LocalDb
from utils.io import Save, Load, Delete
from utils.log import get_logger


logger = get_logger(__name__)


global_finish_sn_name = "Finish_Sn"


def get_young_diagrams_file_name(s_n: int, is_full_path=False):
    """
    Sn or <top_path>/cgc_results/young_diagrams_info/Sn
    p.s. S7 or <top_path>/cgc_results/young_diagrams_info/S7
    """
    from conf.cgc_config import young_diagrams_file_name_format
    if not isinstance(s_n, int) or not isinstance(is_full_path, bool):
        err_msg = "s_n={} with type={} must be int, is_full_path={} with type={} must be bool".format(
            s_n, type(s_n), is_full_path, type(is_full_path))
        logger.error(err_msg)
        return False, err_msg

    file_name = young_diagrams_file_name_format.format(str(s_n))
    if not is_full_path:
        return True, file_name
    else:
        full_path = os.path.join(top_path, cgc_rst_folder, "young_diagrams_info", file_name)
        return True, full_path


def get_young_diagrams_finish_s_n_name(is_full_path=False):
    if not isinstance(is_full_path, bool):
        err_msg = "is_full_path={} with type={} must be bool".format(is_full_path, type(is_full_path))
        logger.error(err_msg)
        return False, err_msg

    if not is_full_path:
        return True, global_finish_sn_name
    else:
        full_path = os.path.join(top_path, cgc_rst_folder, "young_diagrams_info", global_finish_sn_name)
        return True, full_path


class CGCLocalDb(LocalDb):
    """
    这是一个为了CGC特化的父类
    主要修改了init里的参数，增加了txt接口用于展示
    也需要被同目录下的typing继承，才能执行完整功能
    注意：
    db的交互依然是走pkl
    txt展示，应该使用额外的函数，避免逻辑混乱
    """

    def __init__(self):
        """定义connect路径，必要参数等等"""
        super(CGCLocalDb, self).__init__()
        self.db_folder = os.path.join(top_path, cgc_rst_folder, "{}", "{}")  # <top_path>/cgc_results/{}/{}
        self.design_table_type.update({
            "file_name": str,
            "flags": None
        })  # db会自动添加create_time:str和last_write_time:str两项
        self.table_type = None  # 标志数据文件夹意义的静态名称如"young_diagram_info"
        self.map_id = "file_name"
        self.s_n = 0
        self.txt_limit = 0  # 来自config，用来判断txt格式存到S几

    def _get_file_path_without_type_by_condition(self, condition: dict):
        # 只管拼，不管检查，使用者自行检查
        # <top_path>/cgc_results/xxx_info/file_name
        # file_name中是可以包含多层次级目录的，但是不包含文件扩展名
        return super(CGCLocalDb, self)._get_file_path_by_condition(condition, file_type="")

    def _init_cgc_static_db_folder(self):
        """供类继承者创建静态目录，到xxx_info级别"""
        super(CGCLocalDb, self)._init_db_folder()

    def query_by_file_name(self, file_name):
        return self.query({self.map_id: file_name})

    def delete_by_file_name(self, file_name):
        return self.delete({self.map_id: file_name})

    def update_by_file_name(self, file_name, partial_table):
        return self.update({self.map_id: file_name}, partial_table)

    def _query_txt(self, condition: dict, is_auto_limit=True):
        if not isinstance(is_auto_limit, bool):
            err_msg = "is_auto_limit={} must be bool".format(is_auto_limit)
            logger.error(err_msg)
            return False, err_msg
        if is_auto_limit and self.txt_limit != -1 and self.s_n > self.txt_limit:
            # 自动判断 并且判断结果为应该限制
            logger.debug("limit txt with is_auto_limit={}, s_n={}, txt_limit".format(
                is_auto_limit, self.s_n, self.txt_limit))
            return True, False

        condition_copy = copy.deepcopy(condition)
        # query txt 不检查数据
        # txt 里面存的，是table的一个键的值
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(condition_copy)
        if not flag:
            return flag, file_path_without_type  # 此时的file_path_without_type是err_msg
        file_path_txt = file_path_without_type + ".txt"  # query目前操作pkl
        if not os.path.exists(file_path_txt):
            # debug太多了 先注释掉
            # err_msg = "file_path_txt={} not existed, will not load, pls check".format(file_path_txt)
            # logger.debug(err_msg)  # query不需要给warning，delete和update可以
            return True, False
        flag, data = Load.load_txt(file_path_txt)
        return flag, data

    def query_txt_by_file_name(self, file_name, is_auto_limit=True):
        return self._query_txt({self.map_id: file_name}, is_auto_limit=is_auto_limit)

    def _delete_txt(self, condition, is_auto_limit=True):
        if not isinstance(is_auto_limit, bool):
            err_msg = "is_auto_limit={} must be bool".format(is_auto_limit)
            logger.error(err_msg)
            return False, err_msg
        if is_auto_limit and self.txt_limit != -1 and self.s_n > self.txt_limit:
            # 自动判断 并且判断结果为应该限制
            logger.debug("limit txt with is_auto_limit={}, s_n={}, txt_limit".format(
                is_auto_limit, self.s_n, self.txt_limit))
            return True, False

        condition_copy = copy.deepcopy(condition)
        flag, data = self._query_txt(condition_copy, is_auto_limit=is_auto_limit)  # 对于condition的常规检查，query里也有，所以本函数中可以免去
        if not flag:
            return flag, data  # 此时data是errmsg
        if flag and (data is False):
            err_msg = "delete condition_copy={} not existed, return True, False".format(condition_copy)
            logger.warning(err_msg)
            logger.debug("it maybe a wrong event with unexpected file_path,"
                         "or maybe a right event with only find no result")
            return True, False
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(condition_copy)
        if not flag:
            return flag, file_path_without_type

        flag, msg = Delete.delete_txt(file_path_without_type)
        if not flag:
            return flag, msg
        return True, True

    def delete_txt_by_file_name(self, file_name, is_auto_limit=True):
        return self._delete_txt({self.map_id: file_name}, is_auto_limit=is_auto_limit)

    def insert_txt(self, table, is_auto_limit=True, point_key="data"):
        if not isinstance(is_auto_limit, bool) or point_key not in self.design_table_type.keys():
            err_msg = "is_auto_limit={} must be bool, and point_key={} must in design_table_type.keys={}".format(
                is_auto_limit, point_key, list(self.design_table_type.keys()))
            logger.error(err_msg)
            return False, err_msg
        if is_auto_limit and self.txt_limit != -1 and self.s_n > self.txt_limit:
            # 自动判断 并且判断结果为应该限制
            logger.debug("limit txt with is_auto_limit={}, s_n={}, txt_limit".format(
                is_auto_limit, self.s_n, self.txt_limit))
            return True, False

        # insert要根据design_table_type检查table
        table_copy = copy.deepcopy(table)
        flag, query_data = self._query_txt(table_copy, is_auto_limit=is_auto_limit)
        if not flag:
            return flag, query_data  # 此时table是errmsg
        if flag and query_data:
            err_msg = "query_data={} has existed, cannot insert table_copy={}".format(query_data, table_copy)
            logger.error(err_msg)
            return False, err_msg

        flag, msg = self._check_table_key_and_value_type_in_design(table_copy, appear_time_return_false=False)
        if not flag:
            return flag, msg
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(table_copy)
        if not flag:
            return flag, file_path_without_type  # 此时的file_path_without_type是err_msg

        # 取data保存
        flag, msg = Save.save_txt(table_copy.get(point_key), file_path_without_type)
        return flag, msg

    def _update_txt(self, condition, partial_table, is_auto_limit=True, point_key="data"):
        """
        只有partial_table的data中了，才更新
        这里update有三种可能结果：
        1，合法：a)找到一个结果并更新；b)找到一个结果，但无需更新。返回True，True
        2，合法：未找到结果，没有更新，但也没错误。返回True，False
        3，非法：有报错。返回False，err_msg
        """
        if not isinstance(is_auto_limit, bool) or point_key not in self.design_table_type.keys():
            err_msg = "is_auto_limit={} must be bool, and point_key={} must in design_table_type.keys={}".format(
                is_auto_limit, point_key, list(self.design_table_type.keys()))
            logger.error(err_msg)
            return False, err_msg
        if is_auto_limit and self.txt_limit != -1 and self.s_n > self.txt_limit:
            # 自动判断 并且判断结果为应该限制
            logger.debug("limit txt with is_auto_limit={}, s_n={}, txt_limit".format(
                is_auto_limit, self.s_n, self.txt_limit))
            return True, False

        condition_copy = copy.deepcopy(condition)
        partial_table_copy = copy.deepcopy(partial_table)
        flag, data = self._query_txt(condition_copy, is_auto_limit=is_auto_limit)  # 对于condition的常规检查，query里也有，所以本函数中可以免去
        if not flag:
            return flag, data  # 此时data是errmsg
        # 检查是否有真的data
        if flag and (data is False):
            err_msg = "update condition_copy={} not existed, return True, False".format(condition_copy)
            logger.warning(err_msg)
            logger.debug("it maybe a wrong event with unexpected file_path,"
                         "or maybe a right event with only find no result")
            return True, False
        # 检查partial_table和data的一致性，如果一致，不需要真实触发更新动作了
        if str(partial_table_copy.get(point_key)) == data:
            logger.debug("txt_data == data in partial_table_copy={}, will return True, True".format(partial_table_copy))
            return True, True
        new_table = copy.deepcopy(condition_copy)
        new_table.update(partial_table_copy)
        # 剩下的就是真正需要更新的了
        # TODO 删和存最好使用原子化操作，最简单的就没考虑那么多，也没加回档，只能前面尽量检查好
        flag, msg = self._delete_txt(condition_copy, is_auto_limit=is_auto_limit)
        if not flag:
            logger.error("DANGEROUS error! a non-atomized operation may happen")
            logger.error("update within delete and insert, but delete fail, insert not do")
            logger.debug("condition_copy={}, partial_table_copy={}".format(condition_copy, partial_table_copy))
            return flag, msg
        flag, msg = self.insert_txt(new_table, is_auto_limit=is_auto_limit, point_key=point_key)
        if not flag:
            logger.error("DANGEROUS error! a non-atomized operation may happen")
            logger.error("update within delete and insert, but delete done, insert fail")
            logger.debug("condition_copy={}, partial_table_copy={}, new_table={}, point_key={}".format(
                condition_copy, partial_table_copy, new_table, point_key))
            return flag, msg
        return True, True

    def update_txt_by_file_name(self, file_name, partial_table, is_auto_limit=True, point_key="data"):
        return self._update_txt({self.map_id: file_name}, partial_table,
                                is_auto_limit=is_auto_limit, point_key=point_key)

    '''
    下面的函数都是存在于父类，但是以防误用，需要在子类中屏蔽掉的
    TODO 以后如果发现更标准的屏蔽方法，替换掉它们
    '''

    def delete_by_id(self):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def query_by_id(self):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def update_by_id(self):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg
