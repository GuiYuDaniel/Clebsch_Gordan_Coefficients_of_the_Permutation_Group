# -*- coding:utf8 -*-
"""
pipeline负责构建（不是执行）整体计算图：
（当前支持串行拓扑序调度）

1，构建顶层拓扑序（既检查了DAG性，又是串行执行顺序）

2，创建pipenodes


本工作目前只考虑单线程串行，所以这是个简略版pipeline，
暂不支持：微服务，并行异步计算
"""

import copy
from db.typing import PipeLineInfo
from pipeline.pipenode import PipeNode
from utils.config import Config
from utils.log import get_logger
from utils.utils import new_id
from utils.topo import (calc_dag, calc_topo_order)


logger = get_logger(__name__)


class PipeLine(object):
    """
    负责检查并构建计算图，分配对应的pipe node
    pipeline在新建后，是不随任务改变的
    """

    def __init__(self, workflow_conf=None, ppl_id=None, ppl_name=None):
        if workflow_conf is None:
            # a blank ppl waiting for load
            self.ppl_id = ppl_id if ppl_id and isinstance(ppl_id, str) else None
            self.ppl_name = ppl_name if ppl_name and isinstance(ppl_name, str) else None
            self.config = None
            self.dag_dict = None  # 计算图的有向无环图表达
            self.topo_order_list = None  # 序列化的计算顺序
            self.node_id_dict = None  # {name: ppd_id}
            self.flags = None
        else:
            # new a real ppl
            self.ppl_id = ppl_id if ppl_id and isinstance(ppl_id, str) else new_id()
            self.ppl_name = ppl_name if ppl_name and isinstance(ppl_name, str) else self.ppl_id
            self.config = Config()
            self.dag_dict = calc_dag(workflow_conf)
            self.topo_order_list = calc_topo_order(self.dag_dict)
            self.node_id_dict = self._create_node_dict(workflow_conf)
            self.flags = None
            self.save_to_db()

    def _from_dict(self, ppl_dict):
        if not isinstance(ppl_dict, dict):
            err_msg = "_from_dict need input dict but {} with {}".format(type(ppl_dict), ppl_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        map_dict = {
            "pipeline_id": "ppl_id",
            "pipeline_name": "ppl_name",
            "dag_dict": "dag_dict",
            "topo_order_list": "topo_order_list",
            "config": "config",
            "node_id_dict": "node_id_dict",
            "flags": "flags"
        }
        self_keys = self.__dict__.keys()
        for map_key in map_dict:
            # if ppl_dict.get(map_key) and map_dict.get(map_key) in self_keys:
            if map_key in ppl_dict and map_dict.get(map_key) in self_keys:
                self.__dict__[map_dict.get(map_key)] = ppl_dict.get(map_key)

    def load_by_id(self, ppl_id):
        """成功则返回通过id装载起的类实例，失败则返回空实例"""
        if not isinstance(ppl_id, str):
            err_msg = "the ppl_id={} must be str, will return blank class".format(ppl_id)
            logger.error(err_msg)
            # 返回未经改动的实例
            return self
        flag, ppl_dict = PipeLineInfo().query_by_id(ppl_id)
        if not flag:
            err_msg = "query pipeline by ppl_id={} meet error with {}, will return blank class".format(
                ppl_id, ppl_dict)  # 此时ppl_dict是msg
            logger.error(err_msg)
            return self
        if ppl_dict is False:
            err_msg = "no pipeline in db by ppl_id={}, will return blank class".format(ppl_id)
            logger.warning(err_msg)
            return self
        self._from_dict(ppl_dict)
        return self

    def _to_dict(self):
        ppl_dict = {
            "pipeline_id": self.ppl_id,
            "pipeline_name": self.ppl_name,
            "dag_dict": self.dag_dict,
            "topo_order_list": self.topo_order_list,
            "config": self.config,
            "node_id_dict": self.node_id_dict,
            "flags": self.flags
        }
        return ppl_dict

    def save_to_db(self):  # 不需要写try，因为db交互层已经写好了保护和log，只需要返回结果即可
        ppl_dict = self._to_dict()
        flag, msg = PipeLineInfo().insert(ppl_dict)
        if not flag:
            err_msg = "save ppl_dict={} to db meet error".format(ppl_dict)
            logger.error(err_msg)
        return flag, msg

    def _create_node_dict(self, workflow_conf):
        """
        这个函数一是补完node信息，并存入db；
        另一个是返回{node_name: node_id}的字典存在于pipeline中作为索引
        """
        # 补完prep_nodes添加进workflow_conf
        padding_workflow_conf = copy.deepcopy(workflow_conf)
        for single_node_id_dict in padding_workflow_conf:
            node_name = single_node_id_dict.get("name")
            single_node_id_dict["prep_nodes"] = self.dag_dict.get(node_name).get("prep_nodes")
        # check
        node_name_list = self.topo_order_list
        if len(node_name_list) != len(padding_workflow_conf):
            err_msg = "node_name_list={} and padding_workflow_conf={} must have same len".format(
                node_name_list, padding_workflow_conf)
            logger.error(err_msg)
            raise Exception(err_msg)
        node_id_dict = {}
        # create nodes
        for single_dict in padding_workflow_conf:  # 结果的字典也是无序化的，依workflow字典初始化node没问题
            if single_dict.get("name") not in node_name_list:
                err_msg = "node_name={} not in topo_order_list={}, conf_dict={}".format(
                    single_dict.get("name"), node_name_list, single_dict)
                logger.error(err_msg)
                raise Exception(err_msg)
            try:
                ppn = PipeNode(conf=single_dict)
                node_id_dict[single_dict.get("name")] = ppn.ppn_id
            except Exception as e:
                err_msg = "create pipenode with conf={} meet error".format(single_dict)
                logger.error(err_msg)
                raise Exception(e)
        if len(node_name_list) != len(node_id_dict):
            err_msg = "node_name_list={} and node_id_dict={} must have same len".format(node_name_list, node_id_dict)
            logger.error(err_msg)
            logger.debug("workflow_conf={}".format(workflow_conf))
            raise Exception(err_msg)
        return node_id_dict
