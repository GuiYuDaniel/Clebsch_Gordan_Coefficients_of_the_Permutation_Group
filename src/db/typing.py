# -*- coding:utf8 -*-
"""
传入的key和类型，写在db_api外，当作一个小conf一起传入，db_api根据传入自行判断接受与否
db会自动创建两个时间键值：create_time，last_write_time
"""


from db.local_db import LocalDb as BaseDb


class PipeTaskInfo(BaseDb):

    def __init__(self):
        super(PipeTaskInfo, self).__init__()
        self.table_type = "pipetask_info"
        self.design_table_type.update({  # db会自动添加create_time:str和last_write_time:str两项
            "pipetask_id": str,

            "pipeline_id": str,
            "finish_node_list": list,
            "first_input_args": None,
            "first_input_kwargs": None,

            "pipetask_status": str,
            "flags": None
        })
        self.map_id = "pipetask_id"
        self._init_db_folder()


class PipeLineInfo(BaseDb):

    def __init__(self):
        super(PipeLineInfo, self).__init__()
        self.table_type = "pipeline_info"
        self.design_table_type.update({  # db会自动添加create_time:str和last_write_time:str两项
            "pipeline_id": str,
            "pipeline_name": str,

            "dag_dict": dict,
            "topo_order_list": list,
            "config": None,
            "node_id_dict": dict,
            "flags": None
        })
        self.map_id = "pipeline_id"
        self._init_db_folder()


class PipeNodeInfo(BaseDb):

    def __init__(self):
        super(PipeNodeInfo, self).__init__()
        self.table_type = "pipenode_info"
        self.design_table_type.update({  # db会自动添加create_time:str和last_write_time:str两项
            "pipenode_id": str,
            "pipenode_name": str,

            "func_des": list,
            "func_str": str,
            "type": str,
            "inputs": list,
            "outputs": list,
            "next_nodes": list,
            "prep_nodes": list,
            "outputs_r": dict,
            "flags": None
        })
        self.map_id = "pipenode_id"
        self._init_db_folder()
