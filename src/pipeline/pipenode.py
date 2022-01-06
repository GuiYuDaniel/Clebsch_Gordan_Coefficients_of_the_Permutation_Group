# -*- coding:utf8 -*-
"""
pipenode负责维护单个计算节点的所有信息，
被希望于当需要使用哪个节点明确时，访问对应ppn就可以了：

1，构建节点信息

2，更新节点信息

3，检查节点内函数实例化
"""

from db.typing import PipeNodeInfo
from utils.log import get_logger
from utils.utils import new_id


logger = get_logger(__name__)


class PipeNode(object):
    """
    pipenode负责记录节点信息，由pipeline负责初始化。是计算到具体节点时，供pipetask访问的对象
    pipeline在新建后，随任务改变的量有：outputs_r
    TODO 暂不考虑动态装载func，且默认所有func都是静态存在于项目中，存在于已有PYTHONPATH里
    """

    def __init__(self, conf=None, ppn_id=None):
        self._ppn_id = ppn_id if ppn_id and isinstance(ppn_id, str) else None
        self._ppn_name = None
        self._func_des = None
        self._func_str = None
        self._type = None
        self._inputs = None
        self._outputs = None
        # self._extra_args = None  # todo
        # self._extra_kwargs = None  # todo
        self._next_nodes = None
        self._prep_nodes = None
        self._outputs_r = None  # 用_r区分实际变量与描述性变量
        self._flags = None
        if conf:
            self._create_with_conf(conf=conf, ppn_id=ppn_id)
            self.check_node()
            _, _ = self.save_to_db()

    '''这个类所有参数都不应该被随意更改'''

    @property
    def ppn_id(self):
        return self._ppn_id

    @ppn_id.setter
    def ppn_id(self, ppn_id):
        if ppn_id is not None and not isinstance(ppn_id, str):
            err_msg = "ppn_id={} must be None or str".format(ppn_id)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._ppn_id = ppn_id

    @property
    def ppn_name(self):
        return self._ppn_name

    @ppn_name.setter
    def ppn_name(self, ppn_name):
        if ppn_name is not None and not isinstance(ppn_name, str):
            err_msg = "name={} must be None or str".format(ppn_name)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._ppn_name = ppn_name

    @property
    def func_des(self):
        return self._func_des

    @func_des.setter
    def func_des(self, func):  # 在ppn里存func描述，便于落盘需要
        """
        func形如：['None.xx', 'yy', 'zz'] 或 ['xx', 'yy', 'zz']，
        使得后面代码可以造出form xx import yy as zz
        """
        if func is not None \
                and (not isinstance(func, list) or len(func) != 3 or
                     any(not isinstance(i, str) for i in func)):
            err_msg = "func={} must None or a len 3 list, like ['None.xx', 'yy', 'zz'] or ['xx', 'yy', 'zz']" \
                      "in order translate to: form xx import yy as zz".format(func)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._func_des = func

    @property
    def func_str(self):
        return self._func_str

    def _set_func_str(self):
        # 暂时保留None的设定
        if self.func_des[0].startswith("None"):
            tmp_func_des_0 = self.func_des[0].replace("None.", "")
        else:
            tmp_func_des_0 = self.func_des[0]
        if self.func_des[2]:
            func_str = "from {} import {} as {}".format(tmp_func_des_0, self.func_des[1], self.func_des[2])
        else:
            func_str = "from {} import {}".format(tmp_func_des_0, self.func_des[1])
        self._func_str = func_str

    def check_node(self):
        self._check_func()

    def _check_func(self):
        try:
            exec(self._func_str)
            if self.func_des[2]:
                func_r = eval(self.func_des[2])
            else:
                func_r = eval(self.func_des[1])
            del func_r
            logger.debug("func_r={} checked by import and del with node name={}".format(self._func_str, self._ppn_name))
        except Exception as e:
            logger.error(e)
            err_msg = "func_r={} cannot imported with node_name={}".format(self._func_str, self._ppn_name)
            logger.error(err_msg)
            raise Exception(err_msg)

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type_):
        if type_ is not None and type_ != "cold":
            err_msg = "type={} now only support ['cold']".format(type_)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._type = "cold"

    @property
    def inputs(self):
        return self._inputs

    @inputs.setter
    def inputs(self, inputs):
        if inputs is not None \
                and (not isinstance(inputs, list) or not inputs or any(not isinstance(i, str) for i in inputs)):
            err_msg = "inputs={} must be a real str list with node name={}".format(inputs, self._ppn_name)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._inputs = inputs

    @property
    def outputs(self):
        return self._outputs

    @outputs.setter
    def outputs(self, outputs):
        if outputs is not None \
                and (not isinstance(outputs, list) or not outputs or any(not isinstance(i, str) for i in outputs)):
            err_msg = "outputs={} must be a real list with node name={}".format(outputs, self._ppn_name)
            logger.error(err_msg)
            raise Exception(err_msg)
        self._outputs = outputs

    # @property
    # def extra_args(self):
    #     return self._extra_args
    #
    # @property
    # def extra_kwargs(self):
    #     return self._extra_kwargs

    @property
    def next_nodes(self):
        return self._next_nodes

    @next_nodes.setter
    def next_nodes(self, next_nodes):
        # next_nodes要求是列表，元素是str
        if next_nodes is not None \
                and (not isinstance(next_nodes, list) or any(not isinstance(i, str) for i in next_nodes)):
            err_msg = "next_nodes={} can only be None or str list"
            logger.error(err_msg)
            raise Exception(err_msg)
        self._next_nodes = next_nodes

    @property
    def prep_nodes(self):
        return self._prep_nodes

    @prep_nodes.setter
    def prep_nodes(self, prep_nodes):
        # prep_nodes要求是列表，元素是str
        if prep_nodes is not None \
                and (not isinstance(prep_nodes, list) or any(not isinstance(i, str) for i in prep_nodes)):
            err_msg = "prep_nodes={} can only be None or str list"
            logger.error(err_msg)
            raise Exception(err_msg)
        self._prep_nodes = prep_nodes

    @property
    def flags(self):
        return self._flags

    @flags.setter
    def flags(self, flags):
        # 给未来临时塞东西留个地方，目前没啥要求
        self._flags = flags

    @property
    def outputs_r(self):
        return self._outputs_r

    @outputs_r.setter
    def outputs_r(self, outputs_r):
        # 参数命名方式和存取枚举:
        # 命名方式  例子               存                     取
        # 三段形式  p1:::flag:f_name  {"p1:::flag": value}  {"f_name": value}
        # 省略尾段  p1:::flag         {"p1:::flag": value}  {"flag": value}
        # 省略前缀  flag:f_name       {"flag": value}       {"f_name": value}
        # 一段形式  flag              {"flag": value}       {"flag": value}
        # 本函数返回原则：
        # 保存参数名：前缀+管道参数 > 管道参数
        # 使用参数名：函数参数 > 管道参数
        # 详见 CGC_of_Sn/src/pipeline/pipetask.py:_analysis_param_name
        if outputs_r is not None and not isinstance(outputs_r, dict):
            err_msg = "outputs_r={} must be None or dict".format(outputs_r)
            logger.error(outputs_r)
            raise Exception(err_msg)
        self._outputs_r = outputs_r

    def _from_dict(self, ppn_dict):
        if not isinstance(ppn_dict, dict):
            err_msg = "_from_dict need input dict but {} with {}".format(type(ppn_dict), ppn_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        # todo 用__dict__装载，会导致绕过setter检查，在找到解决办法之前先用if
        # map_dict = {
        #     "pipenode_id": "ppn_id",
        #     "pipenode_name": "ppn_name",
        #
        #     "func_des": "func_des",
        #     "func_str": "func_str",
        #     "type": "type",
        #     "inputs": "inputs",
        #     "outputs": "outputs",
        #     "next_nodes": "next_nodes",
        #     "prep_nodes": "prep_nodes",
        #     "flags": "flag",
        #     "outputs_r": "outputs_r"
        # }
        # self_keys = self.__dict__.keys()
        # for map_key in map_dict:
        #     if ppn_dict.get(map_key) and map_dict.get(map_key) in self_keys:
        #         self.__dict__[map_dict.get(map_key)] = ppn_dict.get(map_key)
        ppn_keys = ppn_dict.keys()
        self.ppn_id = ppn_dict.get("pipenode_id") if "pipenode_id" in ppn_keys else None
        self.ppn_name = ppn_dict.get("pipenode_name") if "pipenode_name" in ppn_keys else None
        self.func_des = ppn_dict.get("func_des") if "func_des" in ppn_keys else None
        self._func_str = ppn_dict.get("func_str") if "func_str" in ppn_keys else None
        self.type = ppn_dict.get("type") if "type" in ppn_keys else None
        self.inputs = ppn_dict.get("inputs") if "inputs" in ppn_keys else None
        self.outputs = ppn_dict.get("outputs") if "outputs" in ppn_keys else None
        self.next_nodes = ppn_dict.get("next_nodes") if "next_nodes" in ppn_keys else None
        self.prep_nodes = ppn_dict.get("prep_nodes") if "prep_nodes" in ppn_keys else None
        self.outputs_r = ppn_dict.get("outputs_r") if "outputs_r" in ppn_keys else None
        self.flags = ppn_dict.get("flags") if "flags" in ppn_keys else None

    def load_by_id(self, ppn_id):
        """成功则返回通过id装载起的类实例，失败则返回空实例"""
        if not isinstance(ppn_id, str):
            err_msg = "the ppn_id={} must be str, will return blank class".format(ppn_id)
            logger.error(err_msg)
            # 返回未经改动的实例
            return self
        flag, ppn_dict = PipeNodeInfo().query_by_id(ppn_id)
        if not flag:
            err_msg = "query pipeline by ppn_id={} meet error with {}, will return blank class".format(
                ppn_id, ppn_dict)  # 此时ppn_dict是msg
            logger.error(err_msg)
            return self
        if ppn_dict is False:
            err_msg = "no pipeline in db by ppn_id={}, will return blank class".format(ppn_id)
            logger.warning(err_msg)
            return self
        self._from_dict(ppn_dict)
        return self

    def _to_dict(self):
        ppn_dict = {
            "pipenode_id": self.ppn_id,
            "pipenode_name": self.ppn_name,

            "func_des": self.func_des,
            "func_str": self.func_str,
            "type": self.type,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "next_nodes": self.next_nodes,
            "prep_nodes": self.prep_nodes,
            "flags": self.flags,
            "outputs_r": self.outputs_r
        }
        return ppn_dict

    def save_to_db(self):  # 不需要写try，因为db交互层已经写好了保护和log，只需要返回结果即可
        ppn_dict = self._to_dict()
        flag, msg = PipeNodeInfo().insert(ppn_dict)
        if not flag:
            err_msg = "save ppn_dict={} to db meet error".format(ppn_dict)
            logger.error(err_msg)
        return flag, msg

    def _create_with_conf(self, conf=None, ppn_id=None):
        if not isinstance(conf, dict):
            err_msg = "create pipenode need real dict conf={}".format(conf)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.ppn_id = ppn_id if ppn_id and isinstance(ppn_id, str) else new_id()
        self.ppn_name = conf.get("name")
        self.func_des = conf.get("func")
        self._set_func_str()  # self.func_str
        self.type = conf.get("type")
        self.inputs = conf.get("inputs")
        self.outputs = conf.get("outputs")
        # self.extra_args = None  # todo
        # self.extra_kwargs = None  # todo
        self.next_nodes = conf.get("next_nodes")
        self.prep_nodes = conf.get("prep_nodes")
        self.outputs_r = {}
        self.flags = conf.get("flags")
