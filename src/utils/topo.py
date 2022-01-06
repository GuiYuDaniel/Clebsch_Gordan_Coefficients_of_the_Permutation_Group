# -*- coding:utf8 -*-
"""
用于DAG到拓扑序
目前仅供pipeline构建计算图调用，
如果未来被很多地方调用，单独拿出来
"""

import copy
from utils.log import get_logger


logger = get_logger(__name__)


def calc_dag(workflow_conf):
    """
    将workflow config先简化为DAG
    其中DAG形式为{name:{"next_nodes": [], "prep_nodes": []}, ...}
    这个DAG足够我们下一步用来计算topo序了
    :param workflow_conf:
    :return: dict
    """
    dag_dict = {}
    # 检查conf并简化出name和对应的next_nodes
    for single_dict in workflow_conf:
        node_name = single_dict.get("name")
        if not node_name:
            err_msg = "dict in {} must have key name".format(single_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        next_nodes = single_dict.get("next_nodes")
        if not isinstance(next_nodes, list):
            err_msg = "next_nodes in {} must be a list".format(single_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        unlaw_nodes = [i for i in next_nodes if not isinstance(i, str)]
        if unlaw_nodes:
            err_msg = "these next nodes={} in {} must str".format(unlaw_nodes, single_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        dag_dict[node_name] = {"next_nodes": next_nodes, "prep_nodes": []}
    # 根据next_nodes回填出prep_nodes
    node_name_list = list(dag_dict.keys())
    for node_name in dag_dict:
        next_nodes = dag_dict.get(node_name).get("next_nodes")
        for single_next_node in next_nodes:
            if single_next_node not in node_name_list:
                err_msg = "next node={} defined in next_nodes={} do not in all node name list={}".format(
                    single_next_node, next_nodes, node_name_list)
                logger.error(err_msg)
                raise Exception(err_msg)
            dag_dict[single_next_node]["prep_nodes"].append(node_name)
    # 注意，虽然dag有了，但dict本身是无序的
    return dag_dict


def _check_dag_dict(dag_dict):
    if not isinstance(dag_dict, dict):
        err_msg = "dag_dict={} must a dict".format(dag_dict)
        logger.error(err_msg)
        return False, err_msg
    unlaw_list = [i for i in dag_dict.values() if
                  not isinstance(i, dict) or [j for j in i.keys() if j not in ["next_nodes", "prep_nodes"]]]
    if unlaw_list:
        err_msg = "unlaw items={} must a dict and only with keys next_nodes and prep_nodes".format(unlaw_list)
        logger.error(err_msg)
        return False, err_msg
    return True, ""


def calc_topo_order(dag_dict):
    """
    由DAG计算topo序，返回topo序list
    order序算法：
    1，拿出所有没有输入的节点追加入order
    2，在它们的输出中，减去拿走的输入
    3，循环，直到不剩或者剩下环
    :param dag_dict: such as {"f1": {"next_nodes": ["f2"], "prep_nodes": []}, ...}
    :return: topo序（list）
    """
    flag, msg = _check_dag_dict(dag_dict)
    if not flag:
        raise msg
    topo_order_list = []
    dag_dict = copy.deepcopy(dag_dict)  # 函数式函数不应随意修改调用者的变量，所以使用deepcopy复制一个副本
    while dag_dict:  # 这个循环要考虑dag_dict不是合法DAG的情况
        no_input_nodes = [i for i in dag_dict if dag_dict.get(i).get("prep_nodes") == []]
        if not no_input_nodes:
            err_msg = "topo_order_list={}, and no null list in next dag_dict={}, pls check".format(
                topo_order_list, dag_dict)
            logger.error(err_msg)
            raise Exception(err_msg)
        topo_order_list += no_input_nodes
        for i in no_input_nodes:
            next_nodes = dag_dict.get(i).get("next_nodes", [])
            for j in next_nodes:
                if not dag_dict.get(j):
                    err_msg = "DAG unlaw:\n" \
                              "{} next nodes is {},\n" \
                              "but node {} not in dag={}".format(i, next_nodes, j, dag_dict)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                if not dag_dict.get(j).get("prep_nodes"):
                    err_msg = "DAG unlaw:\n" \
                              "{} next nodes is {},\n" \
                              "but node {} in dag={} not have key prep_nodes".format(i, next_nodes, j, dag_dict)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                if i not in dag_dict.get(j).get("prep_nodes"):
                    err_msg = "DAG unlaw:\n" \
                              "{} next nodes is {},\n" \
                              "but node {} in dag={} with prep_nodes={} not have {}".format(
                                  i, next_nodes, j, dag_dict, dag_dict.get(j).get("prep_nodes"), i)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                dag_dict[j]["prep_nodes"].remove(i)
            del dag_dict[i]
    return topo_order_list
