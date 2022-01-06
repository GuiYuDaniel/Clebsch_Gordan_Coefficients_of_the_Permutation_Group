# -*- coding:utf8 -*-
"""
这不是一个严格意义上的数据库，只是使用存取pickle，统一接口形式，满足最简单的增删改查需要
"""


import copy
import os
import time
from utils.config import Config
from utils.io import Save, Load, Delete
from utils.log import get_logger


logger = get_logger(__name__)


class LocalDb(object):
    """
    这是一个需要被同目录下typing继承，才能执行功能的类
    这个类负责对接特定的具体db
    typing负责定义table，供模块外部调用
    TODO: 考虑到功能的区分，未来引入其他db时，检查逻辑也可以摘出来写到其他层
    """

    def __init__(self):
        """定义connect路径，必要参数等等"""
        self.db_folder = os.path.join(Config().top_path, "results", "{}", "{}")
        self.table_type = None
        self.design_table_type = {"create_time": str,
                                  "last_write_time": str}
        self.map_id = None

    def _get_file_path_by_condition(self, condition: dict, file_type=".pkl"):
        # 只管拼，不管检查，使用者自行检查
        if any(i is None for i in [self.table_type, self.map_id, condition.get(self.map_id)]):
            err_msg = "table_type={}, map_id={}, condition[map_id]={} all must not None".format(
                self.table_type, self.map_id, condition.get(self.map_id))
            logger.error(err_msg)
            return False, err_msg
        file_path = self.db_folder.format(self.table_type, condition.get(self.map_id)) + file_type
        return True, file_path

    @staticmethod
    def _check_condition_in_table(condition: dict, table: dict):
        if not isinstance(condition, dict) or not isinstance(table, dict):
            err_msg = "condition={} with type={} and table={} with type={} must be dict".format(
                condition, type(condition), table, type(table))
            logger.error(err_msg)
            return False, err_msg
        for key in condition:  # 就是字典交集等于cond
            if key not in table:
                return False, None
            if condition[key] != table[key]:
                return False, None
        return True, None

    def _check_table_key_and_value_type_in_design(self, table, appear_time_return_false=False):
        """只检查table里有的，不能保证table缺失key"""
        if not isinstance(table, dict):
            err_msg = "table={} must be dict".format(table)
            logger.error(err_msg)
            return False, err_msg
        if not self.design_table_type or not isinstance(self.design_table_type, dict):
            err_msg = "design_table_type={} must be dict".format(self.design_table_type)
            logger.error(err_msg)
            return False, err_msg
        for key in table:
            if key not in self.design_table_type:
                err_msg = "key={} in table={} not in design_table_type with key_list={}".format(
                    key, table, list(self.design_table_type.keys()))
                logger.error(err_msg)
                return False, err_msg
            if self.design_table_type.get(key) is not None \
                    and not isinstance(table[key], self.design_table_type.get(key)):
                err_msg = "key={}, value={} with type={} in table={} not by design_table_type={}".format(
                    key, table[key], type(table[key]), table, self.design_table_type)
                logger.error(err_msg)
                return False, err_msg
            if appear_time_return_false and (key in ["create_time", "last_write_time"]):  # 作用于输入的table检查，检查之后再加
                err_msg = "table={} should not add time by hand, system will auto do it".format(table)
                logger.error(err_msg)
                return False, err_msg
        return True, None

    @staticmethod
    def _add_time(table, add_create_time=False):
        time_now = time.asctime()
        if add_create_time:
            table["create_time"] = time_now
        table["last_write_time"] = time_now

    def _insert(self, table, add_create_time=True):  # 这里加add_create_time是为了实现假的update，使用的权宜之计，对外不可见
        # insert要根据design_table_type检查table
        table_copy = copy.deepcopy(table)
        flag, msg = self._check_table_key_and_value_type_in_design(table_copy, appear_time_return_false=add_create_time)
        if not flag:
            return flag, msg
        # 终于可以存了
        flag, file_path = self._get_file_path_by_condition(table_copy)
        if not flag:
            return flag, file_path  # 此时的file_path是err_msg
        self._add_time(table_copy, add_create_time=add_create_time)
        flag, msg = Save.save_pickle(table_copy, file_path)
        return flag, msg

    def insert(self, table):
        """
        这里insert有二种可能结果：
        1，合法：没有重复，创建成功。返回True，True
        2，非法：有重复，没创建；或者有报错。返回False，err_msg
        """
        return self._insert(table)

    def update_by_id(self, table_id, partial_table):
        return self.update({self.map_id: table_id}, partial_table)

    def update(self, condition, partial_table):
        """
        本update的逻辑是：根据condition找数据，按照partial_table里的内容，只更新对应的，没提到的不变
        另外，对于pickle假扮的db，update实质是删了重新存。真的db会支持真正的update
        这里update有三种可能结果：
        1，合法：找到一个结果并更新。返回True，True
        2，合法：未找到结果，没有更新，但也没错误。返回True，False
        3，非法：有报错。返回False，err_msg
        """
        flag, table = self.query(condition)  # 对于condition的常规检查，query里也有，所有本函数中可以免去
        if not flag:
            return flag, table  # 此时table是errmsg
        if table is False:
            err_msg = "cannot find table by condition={}, will not update".format(condition)
            logger.warning(err_msg)
            logger.debug("it maybe a wrong event with unexpected file_path,"
                         "or maybe a right event with only find no result")
            return True, False
        flag, msg = self._check_table_key_and_value_type_in_design(partial_table, appear_time_return_false=True)
        if not flag:
            return flag, msg
        if partial_table.get(self.map_id) and table.get(self.map_id) != partial_table.get(self.map_id):
            err_msg = "non-support update {}, id={} of table must same as id={} of partial_table".format(
                self.map_id, table.get(self.map_id), partial_table.get(self.map_id))
            logger.error(err_msg)
            return False, err_msg
        table.update(partial_table)
        flag, msg = self._check_table_key_and_value_type_in_design(table, appear_time_return_false=False)  # table里是有时间的
        if not flag:
            return flag, msg
        # 上面的检查要尽量保证不会出现除读盘落盘外的逻辑错误
        # TODO 删和存最好使用原子化操作，最简单的就没考虑那么多，也没加回档，只能前面尽量检查好
        flag, msg = self.delete(condition)
        if not flag:
            logger.warning("DANGEROUS error! a non-atomized operation may happen")
            return flag, msg
        flag, msg = self._insert(table, add_create_time=False)
        if not flag:
            logger.warning("DANGEROUS error! a non-atomized operation may happen")
            return flag, msg
        return True, True

    def delete_by_id(self, table_id):
        return self.delete({self.map_id: table_id})

    def delete(self, condition):
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
        flag, file_path = self._get_file_path_by_condition(condition)
        if not flag:  # query成功了，后面基本不用检查了
            return flag, file_path  # 此时的file_path是err_msg
        flag, msg = Delete.delete_pickle(file_path)
        if not flag:
            return flag, msg
        return True, True

    def query_by_id(self, table_id):
        return self.query({self.map_id: table_id})

    def query(self, condition):
        """
        这里query有三种可能结果：
服务于这里的query是先拿id找符合condition的，所以不应出现多个结果
        1，合法：找到一个结果。返回True，data(dict)
        2，合法：未找到结果。返回True，False
        3，非法：有报错。返回False，err_msg
        """
        flag, file_path = self._get_file_path_by_condition(condition)
        if not flag:
            return flag, file_path  # 此时的file_path是err_msg
        if not os.path.exists(file_path):
            err_msg = "file_path={} not existed, will not load, pls check"
            logger.debug(err_msg)  # query不需要给warning，delete和update可以
            return True, False
        flag, data = Load.load_pickle(file_path)
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
