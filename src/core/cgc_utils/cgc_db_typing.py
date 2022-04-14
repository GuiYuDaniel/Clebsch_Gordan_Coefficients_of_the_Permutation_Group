# -*- coding:utf8 -*-
"""
模仿src/db/typing.py写的
传入的key和类型，写在cgc_db api外，当作一个小conf一起传入，cgc_db api根据传入自行判断接受与否
cgc db会自动创建两个时间键值：create_time，last_write_time
"""


from core.cgc_utils.cgc_local_db import CGCLocalDb
from conf.cgc_config import default_s_n


# class StatisticInfo(CGCLocalDb):
#     """
#     这个db用来标记各项max已计算的Sn
#     file_name就是其他项的table_type
#     data就是{"Sn": s_n}
#     txt_limit = -1 表示全存txt
#     """
#
#     def __init__(self):
#         super(CGCLocalDb, self).__init__()
#         self.table_type = "calculated_tables_info"  # <top_path>/cgc_results/calculated_tables_info/file_name
#         self.map_id = "file_name"
#         self.design_table_type.update({
#             "data": int,
#         })
#         self.txt_limit = calculated_tables_txt_limit
#         self._init_cgc_static_db_folder()

class YamanouchiMatrixInfo(CGCLocalDb):
    """
    这个db用来存Yamanouchi对换矩阵（仅限(ij)和(in)）

    Yamanouchi对换矩阵的落盘格式为：
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/ij(i,j).pkl 或
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/in(i,n).pkl     ->
    {
    "file_name": "Sn/[ν_i]/ij(i,j)",
    "data": matrix_ij 或 matrix_in,           # list(list(float))
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    ij表示临近对换;
    in表示末尾对换;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1]/ij(2, 3).pkl: {(2,3): [[-0.5, 0.8660254037844386], [0.8660254037844386, 0.5]]}
    S3/[2, 1]/in(1, 3).pkl: {(1,3): [[-0.5, -0.8660254037844386], [-0.8660254037844386, 0.5]]}

    另存：
    <CG>/yamanouchi_matrix_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": [],
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n}}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典
    """

    def __init__(self, s_n):
        super(YamanouchiMatrixInfo, self).__init__()
        self.table_type = "yamanouchi_matrix_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": list,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class YoungTableInfo(CGCLocalDb):
    """
    这个db用来存杨盘

    杨盘的落盘格式为：
    <CG>/young_tableaux_info/Sn/[ν_i].pkl  ->
    {
    "file_name": "Sn/[ν_i]",
    "data": {"m_i": young_table}
    "flags": {"speed_time": speed_time
              "total_num": total_num}  # 注意，这里的total_num就是[ν_i]_num里的total_num
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    total_num就是len({"m_i": young_table})，表示杨盘总数;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1].pkl: {"1": [[1, 2], [3]], "2": [[1, 3], [2]]}


    另存：
    <CG>/young_tableaux_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": {},
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n}}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典


    另存：
    <CG>/young_tableaux_info/Sn/[ν_i]_num.pkl   ->
    {
    "file_name": "Sn/[ν_i]_num",
    "data": {"total_num": total_num}
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    total_num就是len({"m_i": young_table})，表示杨盘总数;

    例如：
    S3/[2, 1]_num.pkl: {"total_num": 2}
    """

    def __init__(self, s_n):
        super(YoungTableInfo, self).__init__()
        self.table_type = "young_tableaux_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": dict,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class BranchingLawInfo(CGCLocalDb):
    """
    这个db用来存杨图的分支律

    分支律的落盘格式为：
    <top_path>/cgc_results/branching_laws_info/Sn/[ν_i].pkl  ->
    {
    "file_name": "Sn/[ν_i]",
    "data": {
            "BL_num": bl_num,    # int，分支数
            "rows": rows,        # list(int)，去掉格子的py行号列表
            "cols": cols,        # list(int)，去掉格子的py列号列表
            "before_YD": [ν’]    # list(list(int))，前置杨盘列表
            }
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    speed_time表示计算用时（秒）

    例如：
    S3/[2, 1].pkl: {"BL_num": 2, "rows": [1, 0], "cols": [0, 1], "before_YD": [[2], [1, 1]]}


    另存：
    <top_path>/cgc_results/branching_laws_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": {},
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n}}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典
    """

    def __init__(self, s_n):
        super(BranchingLawInfo, self).__init__()
        self.table_type = "branching_laws_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": dict,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class YoungDiagramInfo(CGCLocalDb):
    """
    这个db用来存杨图的配分

    杨图的落盘格式为：
    <top_path>/cgc_results/young_diagrams_info/Sn.pkl ->
    {
    "file_name": "Sn",
    "data": [[gamma_i], ...]
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [[gamma_i], ...]表示配分，是二维列表;
    speed_time表示计算用时（秒）

    例如：
    S3.pkl: [[3], [2, 1], [1, 1, 1]]
    S4.pkl: [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]


    另存：
    <top_path>/cgc_results/young_diagrams_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": [],
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n}}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典
    """

    def __init__(self, s_n):
        super(YoungDiagramInfo, self).__init__()
        self.table_type = "young_diagrams_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": list,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()
