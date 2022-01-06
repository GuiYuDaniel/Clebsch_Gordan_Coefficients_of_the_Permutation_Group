# -*- coding:utf8 -*-
"""
pipetask负责循环体调度：

起始作业
while循环：
    下一个node直到结束
"""

from db.typing import PipeTaskInfo, PipeLineInfo, PipeNodeInfo
from utils.log import get_logger
from utils.utils import PipeTaskStatus, new_id
from pipeline.pipeline import PipeLine


logger = get_logger(__name__)


class PipeTask(object):
    """
    真正干活儿的类
    TODO 未来如果需要引入微服务的话，过程中就要只抛错不raise了
    """

    def __init__(self, ppl=None, use_name_replace_id=None):
        if ppl is None:
            # a blank ppt waiting for load
            self.ppt_id = None
            self.ppl_id = None
            self.finish_node_list = None
            self.first_input_args = None
            self.first_input_kwargs = None
            self.ppt_status = None
            self.flags = None
        else:
            # new a real ppt
            if not isinstance(ppl, PipeLine):
                err_msg = "the type of ppl={} must be PipeLine".format(type(ppl).__name__)
                logger.error(err_msg)
                raise Exception(err_msg)
            if not isinstance(ppl.ppl_id, str):
                err_msg = "ppl_id={} must be a str like, pls check your input ppl".format(ppl.ppl_id)
                logger.error(err_msg)
                logger.debug("more information to help check, there is the dict of input ppl={}".format(ppl.__dict__))
                raise Exception(err_msg)
            self.ppt_id = use_name_replace_id \
                if use_name_replace_id and isinstance(use_name_replace_id, str) else new_id()
            self.ppl_id = ppl.ppl_id
            self.finish_node_list = []
            self.first_input_args = None
            self.first_input_kwargs = None
            self.ppt_status = PipeTaskStatus.PREPARATION.name
            self.flags = None
            # 建完存一手
            _, _ = self.save_to_db()

    def _from_dict(self, ppt_dict):
        if not isinstance(ppt_dict, dict):
            err_msg = "_from_dict need input dict but {} with {}".format(type(ppt_dict), ppt_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        map_dict = {
            "pipetask_id": "ppt_id",
            "pipeline_id": "ppl_id",
            "finish_node_list": "finish_node_list",
            "first_input_args": "first_input_args",
            "first_input_kwargs": "first_input_kwargs",
            "pipetask_status": "ppt_status",
            "flags": "flags"
        }
        self_keys = self.__dict__.keys()
        for map_key in map_dict:
            # if ppt_dict.get(map_key) and map_dict.get(map_key) in self_keys:
            if map_key in ppt_dict and map_dict.get(map_key) in self_keys:
                self.__dict__[map_dict.get(map_key)] = ppt_dict.get(map_key)

    def load_by_id(self, ppt_id):
        """成功则返回通过id装载起的类实例，失败则返回空实例"""
        if not isinstance(ppt_id, str):
            err_msg = "the ppt_id={} must be str".format(ppt_id)
            logger.error(err_msg)
            # 返回未经改动的实例
            return self
        flag, ppt_dict = PipeTaskInfo().query_by_id(ppt_id)
        if not flag:
            err_msg = "query pipetask by ppt_id={} meet error with {}".format(ppt_id, ppt_dict)  # 此时ppt_dict是msg
            logger.error(err_msg)
            return self
        if ppt_dict is False:
            err_msg = "no pipetask in db by ppt_id={}, will return blank class".format(ppt_id)
            logger.warning(err_msg)
            return self
        self._from_dict(ppt_dict)
        return self

    def _to_dict(self):
        ppt_dict = {
            "pipetask_id": self.ppt_id,

            "pipeline_id": self.ppl_id,
            "finish_node_list": self.finish_node_list,
            "first_input_args": self.first_input_args,
            "first_input_kwargs": self.first_input_kwargs,
            "pipetask_status": self.ppt_status,
            "flags": self.flags
        }
        return ppt_dict

    def save_to_db(self):  # 不需要写try，因为db交互层已经写好了保护和log，只需要返回结果即可
        ppt_dict = self._to_dict()
        flag, msg = PipeTaskInfo().insert(ppt_dict)
        if not flag:
            err_msg = "save ppn_dict={} to db meet error".format(ppt_dict)
            logger.error(err_msg)
        return flag, msg

    def _simple_start(self, *args, restart=False, **kwargs):
        """
        这是一个按列表顺序执行的最简单的调度系统
        业务逻辑的错误，返回FAIL和msg
        pipetask框架的错误，使用raise断掉
        """
        # TODO 添加状态
        # TODO 想清楚执行单元，是回访db获得信息，还是从类获得信息，更合理。是否可能出现不一致
        # 检查ppt和ppl，获得状态，进度
        if not self.ppl_id or not isinstance(self.ppl_id, str) \
                or not self.ppt_id or not isinstance(self.ppt_id, str):
            err_msg = "start a pipetask must has real pipeline_id and pipetask_id, " \
                      "with ppl_id={}, ppt_id={}".format(self.ppl_id, self.ppt_id)
            logger.error(err_msg)
            logger.debug("full ppt_dict for debug {}".format(self._to_dict()))
            raise Exception(err_msg)
        flag, ppl_info = PipeLineInfo().query_by_id(self.ppl_id)
        if not flag or not ppl_info:
            err_msg = "cannot get pipeline from db, with pipeline_id={}, error_msg={}".format(self.ppl_id, ppl_info)
            logger.error(err_msg)
            raise Exception(err_msg)
        topo_order_list = ppl_info.get("topo_order_list")
        node_id_dict = ppl_info.get("node_id_dict")
        if not topo_order_list or not isinstance(topo_order_list, list):
            err_msg = "topo_order_list={} must be a not null list".format(topo_order_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        if not node_id_dict or not isinstance(node_id_dict, dict):
            err_msg = "node_id_dict={} must be a not null list".format(node_id_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        if len(topo_order_list) != len(node_id_dict):
            err_msg = "topo_order_list={} and node_id_dict={} must have same len".format(topo_order_list, node_id_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        if any(not node_id_dict.get(i) or not isinstance(node_id_dict.get(i), str) for i in topo_order_list):
            err_msg = "every name in topo_order_list={} must be key in node_id_dict={}, " \
                      "and the value must be str".format(topo_order_list, node_id_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        if restart:
            if not isinstance(self.finish_node_list, list):
                err_msg = "finish_node_list={} must a list for restart".format(self.finish_node_list)
                logger.error(err_msg)
                raise Exception(err_msg)
            # 开始按照topo_order_list里的顺序派发作业，但finish_node_list里有的要跳过
            logger.info("\n" + "=" * 80 + "\n" + "=" * 80 + "\n" + "=" * 80)
            logger.info("pipetask will restart by topo_order_list={}, and finish_node_list={} will be jump".format(
                topo_order_list, self.finish_node_list))
            logger.info("\n" + "=" * 80)
            self._transfer_and_update_status_to(PipeTaskStatus.DOING.name,
                                                now_status_check=PipeTaskStatus.RESTARTING.name)
        else:
            # 开始按照topo_order_list里的顺序派发作业
            logger.info("\n" + "=" * 80 + "\n" + "=" * 80 + "\n" + "=" * 80)
            logger.info("pipetask will start by topo_order_list={}".format(topo_order_list))
            logger.info("\n" + "=" * 80)
            self._transfer_and_update_status_to(PipeTaskStatus.DOING.name,
                                                now_status_check=PipeTaskStatus.PREPARATION.name)
        # 真正的节点循环开始！
        for now_node_name in topo_order_list:
            if restart:
                if now_node_name in self.finish_node_list:
                    logger.info("\n" + "=" * 40)
                    logger.info("jump finished pipenode={} for restart".format(now_node_name))
                    logger.info("\n" + "=" * 40)
                    logger.info("\n" + "=" * 40)
                    continue
            now_node_id = node_id_dict.get(now_node_name)
            logger.info("\n"+"=" * 40)
            logger.info("pipenode={} start".format(now_node_name))
            logger.info("\n"+"=" * 40)
            try:
                flag, ppn_info = PipeNodeInfo().query_by_id(now_node_id)
                if not flag or not ppn_info:
                    err_msg = "cannot get pipenode from db, with pipenode_id={}, error_msg={}".format(
                        now_node_id, ppn_info)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                func_r = self._get_func_r(ppn_info)
                inputs_r_args, inputs_r_kwargs = self._get_inputs_r(ppn_info, node_id_dict, *args, **kwargs)
                logger.info("calling function={} with args={}, kwargs={} in node={}".format(
                    func_r.__name__, inputs_r_args, inputs_r_kwargs, now_node_name))
                outputs_r = func_r(*inputs_r_args, **inputs_r_kwargs)  # 注意，执行结果有可能很复杂，要尽量测试到所有情况
                logger.info("function={} return {} in node={}".format(func_r.__name__, outputs_r, now_node_name))
                flag, msg = self._update_outputs_to_node(ppn_info, outputs_r)
                if flag is not True or msg is not True:
                    err_msg = "update outputs_r={} to node={} meet fail".format(outputs_r, ppn_info)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                self.finish_node_list.append(now_node_name)
                flag, msg = self._update_pipetask({"finish_node_list": self.finish_node_list})
                if flag is not True or msg is not True:
                    err_msg = "update finish_node_list={} to pipetask={} meet fail".format(
                        self.finish_node_list, self.ppt_id)
                    logger.error(err_msg)
                    raise Exception(err_msg)
            except Exception as e:
                logger.error(e)
                logger.error("pipeline error with node_name={}, function={}".format(now_node_name, func_r.__name__))
                self._transfer_and_update_status_to(PipeTaskStatus.FAIL.name, now_status_check=PipeTaskStatus.DOING.name)
                return self.ppt_status, e
            logger.info("\n"+"=" * 40)
            logger.info("pipenode={} done with outputs={}".format(now_node_name, outputs_r))
            logger.info("pipenode={} succ".format(now_node_name))
            logger.info("now finish_node_list={}".format(self.finish_node_list))
            logger.info("\n"+"=" * 40)
        logger.info("\n" + "=" * 80 + "\n" + "=" * 80 + "\n" + "=" * 80)
        logger.info("All pipeline done with finish_node_list={}".format(self.finish_node_list))
        logger.info("\n" + "=" * 80)
        logger.info("完结撒花！")
        logger.info("～(￣▽￣～) <(￣︶￣)>(～￣▽￣)～")
        self._transfer_and_update_status_to(PipeTaskStatus.SUCCESS.name, now_status_check=PipeTaskStatus.DOING.name)
        return self.ppt_status, None

    def start(self, *args, mode="simple", **kwargs):
        self.first_input_args = args
        self.first_input_kwargs = kwargs
        update_data = {"first_input_args": self.first_input_args,
                       "first_input_kwargs": self.first_input_kwargs}
        flag, msg = self._update_pipetask(update_data)
        if flag is not True or msg is not True:
            # 转译状态失败只抛错不再转移，不然就是无限循环了
            err_msg = "update pipetask_status={} meet fail".format(self.ppt_status)
            logger.error(err_msg)
            raise Exception(err_msg)
        mode_list = ["simple"]
        if mode == "simple":
            return self._simple_start(*args, **kwargs)
        else:
            err_msg = "cannot start a pipetask={} for mode={} not in support mode_list={}".format(
                self.ppt_id, mode, mode_list)
            logger.error(err_msg)
            raise Exception(err_msg)

    def _simple_restart(self):
        """
        1，检查是否可以restart，如果可以，更新状态以及修改db
        2，交给start函数完成后续步骤，start对于完成的节点，采用跳过模式，便于未来加并行代码
        """
        # check
        # check ppt
        if self.ppt_id is None:
            err_msg = "cannot restart a pipetask that pipetask_id is None"
            logger.error(err_msg)
            raise Exception(err_msg)
        other_check_params_list = [self.ppl_id, self.finish_node_list, self.first_input_args, self.first_input_kwargs]
        if any(i is None for i in other_check_params_list):
            err_msg = "this pipetask param list={} must not None for restart pipetask".format(other_check_params_list)
            logger.error(err_msg)
            raise Exception(err_msg)
        restartable_status_list = [PipeTaskStatus.DOING.name, PipeTaskStatus.FAIL.name]
        if self.ppt_status not in restartable_status_list:
            err_msg = "pipetask only support be restart in these status={}, but now status={}".format(
                restartable_status_list, self.ppt_status)
            logger.error(err_msg)
            raise Exception(err_msg)
        # check ppl
        flag, ppl_info = PipeLineInfo().query_by_id(self.ppl_id)
        if not flag or not ppl_info:
            err_msg = "cannot get pipeline from db for restart, with pipeline_id={}, error_msg={}".format(
                self.ppl_id, ppl_info)
            logger.error(err_msg)
            raise Exception(err_msg)
        ppl_check_param_keys_list = ["dag_dict", "topo_order_list", "config", "node_id_dict"]
        if any(ppl_info.get(i) is None for i in ppl_check_param_keys_list):
            err_msg = "this pipeline param list={} must not None for restart pipetask".format(
                [ppl_info.get(i) for i in ppl_check_param_keys_list])
            logger.error(err_msg)
            raise Exception(err_msg)
        # check finish ppn
        node_id_dict = ppl_info.get("node_id_dict")
        for finish_ppn_name in self.finish_node_list:
            finish_ppn_id = node_id_dict.get(finish_ppn_name)
            flag, ppn_info = PipeNodeInfo().query_by_id(finish_ppn_id)
            if not flag or not ppn_info:
                err_msg = "cannot get pipenode from db for restart, with pipenode_id={}, error_msg={}".format(
                    finish_ppn_id, ppn_info)
                logger.error(err_msg)
                raise Exception(err_msg)
            outputs_r = ppn_info.get("outputs_r")
            if not outputs_r or not isinstance(outputs_r, dict):
                err_msg = "finish node={} id={} should have real outputs_r={} but not, cannot restart".format(
                    finish_ppn_name, finish_ppn_id, outputs_r)
                logger.error(err_msg)
                raise Exception(err_msg)
        # change db
        # 目前，只需修改状态就可以了，未完成节点没有东西需要清理
        self._transfer_and_update_status_to(PipeTaskStatus.RESTARTING.name)
        # start
        args = self.first_input_args
        kwargs = self.first_input_kwargs
        return self._simple_start(*args, restart=True, **kwargs)

    def restart(self, mode="simple"):
        """
        1，检查是否可以restart，如果可以，更新状态以及修改db
        2，交给start函数完成后续步骤，start对于完成的节点，采用跳过模式，可便于未来加并行代码
                                                      ^^^^^^^
        """
        mode_list = ["simple"]
        if mode == "simple":
            return self._simple_restart()
        else:
            err_msg = "cannot restart a pipetask={} for mode={} not in support mode_list={}".format(
                self.ppt_id, mode, mode_list)
            logger.error(err_msg)
            raise Exception(err_msg)

    def _get_func_r(self, ppn_info):
        # 获取当前node的function
        func_str = ppn_info.get("func_str")
        func_des = ppn_info.get("func_des")
        exec(func_str)  # 把检查放在pipenode里了，在ppnode里没报错的话，这里肯定通过
        if func_des[2]:
            func_r = eval(func_des[2])
        else:
            func_r = eval(func_des[1])
        return func_r

    @staticmethod
    def _analysis_param_name(full_name):
        # 参数命名方式和存取枚举:
        # 命名方式  例子               存                     取
        # 三段形式  p1:::flag:f_name  {"p1:::flag": value}  {"f_name": value}
        # 省略尾段  p1:::flag         {"p1:::flag": value}  {"flag": value}
        # 省略前缀  flag:f_name       {"flag": value}       {"f_name": value}
        # 一段形式  flag              {"flag": value}       {"flag": value}
        # 本函数返回原则：
        # 保存参数名：前缀+管道参数 > 管道参数
        # 使用参数名：函数参数 > 管道参数
        if not full_name or not isinstance(full_name, str):
            err_msg = "param name={} error, it must be str".format(full_name)
            logger.error(err_msg)
            return None, None
        split_name_list = full_name.split(":")
        if len(split_name_list) == 5 and all(split_name_list[i] for i in [0, 3, 4]):
            # 三段式
            return "{}:::{}".format(split_name_list[0], split_name_list[3]), "{}".format(split_name_list[4])
        elif len(split_name_list) == 2 and all(split_name_list[i] for i in [0, 1]):
            # 省略前缀
            return split_name_list[0], split_name_list[1]
        elif len(split_name_list) == 4 and all(split_name_list[i] for i in [0, 3]):
            # 省略函数参数名
            return "{}:::{}".format(split_name_list[0], split_name_list[3]), "{}".format(split_name_list[3])
        elif len(split_name_list) == 1 and all(split_name_list[i] for i in [0]):
            # 一段式
            return split_name_list[0], split_name_list[0]
        else:
            # 命名错误
            err_msg = "param name={} error, it must be write by rule".format(full_name)
            logger.error(err_msg)
            return None, None

    def _get_inputs_r(self, ppn_info, node_id_dict, *args, **kwargs):
        # 获取当前函数的inputs
        # 起始节点返回输入的inputs_r，其他节点从outputs里找，要求都是函数式函数
        node_name = ppn_info.get("pipenode_name")
        inputs_name_list = ppn_info.get("inputs")
        prep_nodes = ppn_info.get("prep_nodes")
        if not prep_nodes:
            logger.debug("node={} is the first node".format(node_name))
            return args, kwargs
        else:
            # 后面不管必要参数还是可选参数，全部按照dict的形式传
            rst = {}
            # 这一步把前置节点的结果拿过来
            all_prep_outputs_r = {}
            for prep_node_name in prep_nodes:
                prep_node_id = node_id_dict.get(prep_node_name)
                flag, prep_ppn_info = PipeNodeInfo().query_by_id(prep_node_id)
                if not flag or not prep_ppn_info:
                    err_msg = "cannot get pipenode from db, with pipenode_id={}, error_msg={}".format(
                        prep_ppn_info, ppn_info)
                    logger.error(err_msg)
                    return None, None
                prep_outputs_r = prep_ppn_info.get("outputs_r")
                if prep_outputs_r:
                    all_prep_outputs_r.update(prep_outputs_r)
                else:
                    err_msg = "Cannot get prep_outputs_r for now_node={}, prep_node={}, prep_node_id={}".format(
                        node_name, prep_node_name, prep_node_id)
                    logger.error(err_msg)
                    return None, None
            # 这一步开始组装本节点需要的入参
            for input_name in inputs_name_list:
                save_name, use_name = self._analysis_param_name(input_name)
                if save_name in all_prep_outputs_r:
                    rst.update({use_name: all_prep_outputs_r.get(save_name)})
                else:
                    err_msg = "input_name={} not in all_prep_outputs_r={}, for now_node={}, prep_node={}".format(
                        save_name, all_prep_outputs_r, node_name, prep_nodes)
                    logger.error(err_msg)
                    logger.debug("All prep_node_outputs_r will given bellow, and hope this will helpful")
                    for prep_node_name in prep_nodes:
                        prep_node_id = node_id_dict.get(prep_node_name)
                        flag, prep_ppn_info = PipeNodeInfo().query_by_id(prep_node_id)
                        if not flag or not prep_ppn_info:
                            err_msg = "cannot get pipenode from db, with pipenode_id={}, error_msg={}".format(
                                prep_ppn_info, ppn_info)
                            logger.error(err_msg)
                            continue
                        prep_outputs = prep_ppn_info.get("outputs")
                        prep_outputs_r = prep_ppn_info.get("outputs_r")
                        logger.debug("node={}, outputs={}, outputs_r={}".format(
                            prep_node_name, prep_outputs, prep_outputs_r))
                    return None, None
            return tuple([]), rst

    def _update_outputs_to_node(self, ppn_info, outputs_r):
        # 对应单输出和多输出，outputs可能是普通单变量或者tuple
        # 但是，单输出也有tuple的情况，要与多变量做区别
        ppn_name = ppn_info.get("pipenode_name")
        ppn_id = ppn_info.get("pipenode_id")
        outputs_name_list = ppn_info.get("outputs")
        # 因为python支持多return，无法根据函数的定义直接判断结果与conf里定义的是否一致，只能根据结果来判断
        if len(outputs_name_list) == 1:
            # 这种情况直接存
            save_name, use_name = self._analysis_param_name(outputs_name_list[0])
            outputs_r_dict = {save_name: outputs_r}
        else:
            # 这种情况存之前要先检查
            if not isinstance(outputs_r, tuple):
                err_msg = "outputs_r={} must be tuple but it is {} when outputs multi with outputs_name={}, \n" \
                          "node={}".format(outputs_r, type(outputs_r), outputs_name_list, ppn_name)
                logger.error(err_msg)
                raise Exception(err_msg)
            if len(outputs_name_list) != len(outputs_r):
                err_msg = "outputs_r={} must same len with outputs_name={}, \n" \
                          "node={} ".format(outputs_r, outputs_name_list, ppn_name)
                logger.error(err_msg)
                raise Exception(err_msg)
            # TODO 目前的检查下，会漏一种情况，就是定义多输出的函数，实际是单输出tuple，且tuple内变量个数与定义的一致
            outputs_r_dict = {}
            for name, value in zip(outputs_name_list, outputs_r):
                save_name, use_name = self._analysis_param_name(name)
                outputs_r_dict[save_name] = value
        update_dict = {"outputs_r": outputs_r_dict}
        flag, msg = PipeNodeInfo().update_by_id(ppn_id, update_dict)
        return flag, msg

    def _transfer_status_to(self, target_status, now_status_check=None):
        if now_status_check:
            if self.ppt_status is not now_status_check:
                err_msg = "pipetask_status should be {} not {}".format(now_status_check, self.ppt_status)
                logger.error(err_msg)
                raise Exception(err_msg)
        next_status_able_list = eval("PipeTaskStatus.{}.value".format(self.ppt_status))
        if target_status not in next_status_able_list:
            err_msg = "status error, now_status={} can transfer to {}, not {}".format(
                self.ppt_status, next_status_able_list, target_status)
            logger.error(err_msg)
            raise Exception(err_msg)
        self.ppt_status = target_status

    def _transfer_and_update_status_to(self, target_status, now_status_check=None):
        # transfer
        self._transfer_status_to(target_status, now_status_check=now_status_check)
        # update
        flag, msg = self._update_pipetask({"pipetask_status": target_status})
        if flag is not True or msg is not True:
            # 转译状态失败只抛错不再转移，不然就是无限循环了
            if now_status_check:
                err_msg = "update pipetask_status from {} to {} meet fail".format(now_status_check, target_status)
            else:
                err_msg = "update pipetask_status={} meet fail".format(target_status)
            logger.error(err_msg)
            raise Exception(err_msg)
        if now_status_check:
            logger.info("change pipetask_id={} status from {} to {}".format(
                self.ppt_id, now_status_check, target_status))
        else:
            logger.info("change pipetask_id={} status to {}".format(self.ppt_id, target_status))

    def _update_pipetask(self, partial_dict, ppt_id=None):
        if not ppt_id:
            flag, msg = PipeTaskInfo().update_by_id(self.ppt_id, partial_dict)
        else:
            logger.debug("save pipetask with other pointed ppt_id={} in ppt_id={} class".format(ppt_id, self.ppt_id))
            flag, msg = PipeTaskInfo().update_by_id(ppt_id, partial_dict)
        return flag, msg
