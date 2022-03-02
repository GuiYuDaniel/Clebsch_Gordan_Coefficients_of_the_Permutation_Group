# -*- coding:utf8 -*-
"""
模仿src/db/typing.py写的
传入的key和类型，写在cgc_db api外，当作一个小conf一起传入，cgc_db api根据传入自行判断接受与否
cgc db会自动创建两个时间键值：create_time，last_write_time
"""


from core.cgc_utils.cgc_local_db import CGCLocalDb
from conf.cgc_config import calculated_tables_txt_limit
from conf.cgc_config import young_diagrams_txt_limit


class CalculatedTableInfo(CGCLocalDb):
    """
    这个db用来标记各项max已计算的Sn
    file_name就是其他项的table_type
    data就是{"Sn": s_n}
    txt_limit = -1 表示全存txt
    """

    def __init__(self):
        super(CGCLocalDb, self).__init__()
        self.table_type = "calculated_tables_info"  # <top_path>/cgc_results/calculated_tables_info/file_name
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": int,
        })
        self.txt_limit = calculated_tables_txt_limit
        self._init_cgc_static_db_folder()


class YoungDiagramInfo(CGCLocalDb):

    def __init__(self, s_n):
        super(CGCLocalDb, self).__init__()
        self.table_type = "young_diagrams_info"  # <top_path>/cgc_results/young_diagrams_info/file_name
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": list,
        })
        self.s_n = s_n
        self.txt_limit = young_diagrams_txt_limit
        self._init_cgc_static_db_folder()
