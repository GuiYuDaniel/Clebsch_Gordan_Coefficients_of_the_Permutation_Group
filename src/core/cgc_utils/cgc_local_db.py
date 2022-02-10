# -*- coding:utf8 -*-
"""
这不是一个严格意义上的数据库，只是使用存取pickle，统一接口形式，满足最简单的增删改查需要
在project的db基础上，多了txt的保存需求，需要改写一些函数
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


class CGCLocalDb(LocalDb):

    def __init__(self):
        super(CGCLocalDb, self).__init__()
        self.db_folder = os.path.join(top_path, cgc_rst_folder, "{}", "{}")  # <top_path>/cgc_results/{}/{}
        self.design_table_type.update({
            "file_name": str,

            "data": dict,
            "flags": None
        })  # db会自动添加create_time:str和last_write_time:str两项
        self.table_type = None  # 标志数据文件夹意义的静态名称如"young_diagram_info"
        self.map_id = "file_name"
        self.s_n = 0
        self.txt_limit = 0  # 来自config，用来判断txt格式存到S几

    def _get_file_path_without_type_by_condition(self, condition: dict):
        # 只管拼，不管检查，使用者自行检查
        # TODO 这个函数后面还需要改造呢，或者在Path里，加个拼装path取代format
        if any(i is None for i in [self.table_type, self.map_id, condition.get(self.map_id)]):
            err_msg = "table_type={}, map_id={}, condition[map_id]={} all must not None".format(
                self.table_type, self.map_id, condition.get(self.map_id))
            logger.error(err_msg)
            return False, err_msg
        # <top_path>/cgc_results/xxx_info/file_name
        # file_name中是可以包含多层次级目录的，但是不包含文件扩展名
        file_path_without_type = self.db_folder.format(self.table_type, condition.get(self.map_id))
        return True, file_path_without_type

    def _init_cgc_static_db_folder(self):
        """供类继承者创建静态目录，到xxx_info级别"""
        super(CGCLocalDb, self)._init_db_folder()

    # TODO 重写一些函数，要考虑深度不定（好办，用个函数拼），类型可判断（类型默认pkl+txt，有强制接口），失败回收
    def _insert(self, table, add_create_time=True, force_file_type=None):
        # insert要根据design_table_type检查table
        table_copy = copy.deepcopy(table)
        flag, msg = self._check_table_key_and_value_type_in_design(table_copy, appear_time_return_false=add_create_time)
        if not flag:
            return flag, msg
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(table_copy)
        if not flag:
            return flag, file_path_without_type  # 此时的file_path_without_type是err_msg
        self._add_time(table_copy, add_create_time=add_create_time)

        # 终于可以存了
        if not force_file_type:  # default mode
            # 保存pkl
            flag, msg = Save.save_pickle(table_copy, file_path_without_type)
            # 按照config决定是否存个txt
            if flag and self.s_n <= self.txt_limit:
                # txt格式只存table_copy["data"]
                flag, msg = Save.save_txt(table_copy.get("data"), file_path_without_type)
                if not flag:
                    # 需要回收pkl
                    logger.warning("save pkl success but txt not! will delete pkl")
                    flag, msg = Delete.delete_pickle(file_path_without_type)
                    if flag:
                        logger.info("delete {}.pkl success".format(file_path_without_type))
                    else:
                        logger.error(msg)

            return flag, msg

        elif force_file_type == ".pkl":
            # 强制保存为pkl
            flag, msg = Save.save_pickle(table_copy, file_path_without_type)
            return flag, msg

        elif force_file_type == ".txt":
            # 强制保存为txt
            flag, msg = Save.save_txt(table_copy, file_path_without_type)
            return flag, msg

        else:
            # 强制保存的格式不支持
            err_msg = "only support .pkl or .txt type, force_file_type={} not supported".format(force_file_type)
            logger.error(err_msg)
            return False, err_msg

    def insert(self, table, force_file_type=None):
        """
        这里insert有二种可能结果：
        1，合法：没有重复，创建成功。返回True，True
        2，非法：有重复，没创建；或者有报错。返回False，err_msg
        """
        return self._insert(table, force_file_type=force_file_type)

    def delete_by_file_name(self, file_name):
        return self.delete({self.map_id: file_name})

    def delete(self, condition, force_file_type=None):
        """
        这里delete有三种可能结果：
        1，合法：找到一个结果并删除。返回True，True
        2，合法：未找到结果，没有删除。返回True，False
        3，非法：有报错。返回False，err_msg
        """
        flag, table = self.query(condition)  # 对于condition的常规检查，query里也有，所以本函数中可以免去
        if not flag:
            return flag, table  # 此时table是errmsg
        if flag and (table is False):
            err_msg = "file_path={} not existed, cannot delete but return True"
            logger.warning(err_msg)
            logger.debug("it maybe a wrong event with unexpected file_path,"
                         "or maybe a right event with only find no result")
            return True, False
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(condition)
        if not flag:
            return flag, file_path_without_type  # 此时file_path_without_type是errmsg

        # 终于可以删了
        if not force_file_type:  # default mode
            # 删除pkl
            flag, msg = Delete.delete_pickle(file_path_without_type)
            # 按照config决定是否删除txt
            if flag and self.s_n <= self.txt_limit:
                flag, msg = Delete.delete_txt(file_path_without_type)
                if not flag:
                    logger.error("delete pkl success but txt not, pls check, with {}".format(file_path_without_type))
                    logger.warning("del pkl but txt, it is dangerous! it may break next inserting! MUST check it!")

            return flag, msg

        elif force_file_type == ".pkl":
            flag, msg = Delete.delete_pickle(file_path_without_type)
            if not flag:
                return flag, msg
            return True, True

        elif force_file_type == ".txt":
            flag, msg = Delete.delete_txt(file_path_without_type)
            if not flag:
                return flag, msg
            return True, True

        else:
            # 强制删除的格式不支持
            err_msg = "only support .pkl or .txt type, force_file_type={} not supported".format(force_file_type)
            logger.error(err_msg)
            return False, err_msg

    def query_by_file_name(self, file_name):
        return self.query({self.map_id: file_name})

    def query(self, condition):
        """
        这里query有三种可能结果：
        服务于这里的query是先拿id找符合condition的，所以不应出现多个结果
        1，合法：找到一个结果。返回True，data(dict)
        2，合法：未找到结果。返回True，False
        3，非法：有报错。返回False，err_msg
        """
        flag, file_path_without_type = self._get_file_path_without_type_by_condition(condition)
        if not flag:
            return flag, file_path_without_type  # 此时的file_path_without_type是err_msg
        file_path_pkl = file_path_without_type + ".pkl"  # query目前操作pkl
        if not os.path.exists(file_path_pkl):
            err_msg = "file_path={} not existed, will not load, pls check"
            logger.debug(err_msg)  # query不需要给warning，delete和update可以
            return True, False
        flag, data = Load.load_pickle(file_path_pkl)
        if not flag:
            return flag, data  # 此时data是errmsg
        # query结果不可信，需要检查
        flag, msg = self._check_table_key_and_value_type_in_design(data, appear_time_return_false=False)
        if not flag:
            return flag, msg
        flag, msg = self._check_condition_in_table(condition, data)
        if not flag:
            if not msg:
                return True, False
            else:
                return False, msg
        return True, data

    '''
    下面的函数都是存在于父类，但是以防误用，需要在子类中屏蔽掉的
    TODO 以后如果发现更标准的屏蔽方法，替换掉它们
    '''
    def _get_file_path_by_condition(self, condition: dict, file_type=".pkl"):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def update_by_id(self, table_id, partial_table):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def update(self, condition, partial_table):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def delete_by_id(self, table_id):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg

    def query_by_id(self, table_id):
        # 相当于删去父类中的这个函数
        err_msg = "function={} cannot be used in class=CGCLocalDb and children class={}".format(
            sys._getframe().f_code.co_name, self.__class__.__name__)
        logger.error(err_msg)
        raise err_msg
