# -*- coding:utf8 -*-
"""
模仿src/db/typing.py写的
传入的key和类型，写在cgc_db api外，当作一个小conf一起传入，cgc_db api根据传入自行判断接受与否
cgc db会自动创建两个时间键值：create_time，last_write_time
"""


import numpy as np
import sympy as sp
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

class CGCInfo(CGCLocalDb):
    """
    这个db用来存CGC

    <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_β_m.pkl  # 无多重性，则为[ν]_m.pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_β_m,  # 无多重性，则为[ν]_m
    "data": cgc_square_dict,
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群；
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    β对应[ν]的多重性；
    m是[ν]的yamanouchi序；
    cgc_square_dict数值是cgc的平方，符号是平方前cgc的符号；
    speed_time表示计算用时（秒）

    例如，
    <CG>/cgc_info/S4/[2, 2]_[3, 1]/[2, 1, 1]_m2.pkl
    {(2, 1): 1/2, (1, 3): 1/4, (2, 2): 1/4}

    <CG>/cgc_info/S5/[3, 1, 1]_[3, 2]/[1, 1, 1]_β1_m4.pkl
    {(1, 4): 1/8, (2, 4): -1/16, (3, 5): 1/16, (4, 1): -1/10, (4, 2): -1/5,
     (4, 4): 3/80, (5, 3): 1/5, (5, 5): -3/80, (6, 3): -1/10, (6, 5): -3/40}
    """
    def __init__(self, s_n):
        super(CGCInfo, self).__init__()
        self.table_type = "cgc_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": (dict, sp.Matrix),
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class EInfo(CGCLocalDb):
    """
    这个db用来存ϵ

    注意，为了方便书写，类名用的是英文字母e的大写，而不是希腊字母ϵ的大写

    <CG>/ϵ_info/Sn/[σ]_[μ]/[ν]_β.pkl  # 无多重性，则为[ν].pkl
    {
    "file_name": Sn/[σ]_[μ]/[ν]_β,  # 无多重性，则为[ν]
    "data": ϵ_dict,
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型；
    β对应[ν]的多重性;

    例如，
    <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 2]_2.pkl
    {"data": {"ϵ0": 1, "ϵ1": 1, "ϵ4": -1, "ϵ14": -1,
                  "ϵ5": 1, "ϵ15": -1, "ϵ6": -1, "ϵ16": 1},
     "flags": {"ϵ0": (1, 1), "ϵ1": (1, 1), "ϵ4": (6, 6), "ϵ14": (6, 6),
               "ϵ5": (6, 4), "ϵ15": (3, 1), "ϵ6": (1, 3), "ϵ16": (4, 6)}}
    <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[2, 1, 1]_2.pkl
    {"data": {"ϵ0": 1, "ϵ1": -1, "ϵ4": -1, "ϵ14": 1,
              "ϵ5": -1, "ϵ15": -1, "ϵ6": 1, "ϵ16": 1},
     "flags": {"ϵ0": (1, 4), "ϵ1": (4, 1), "ϵ4": (6, 3), "ϵ14": (3, 6),
               "ϵ5": (6, 1), "ϵ15": (6, 1), "ϵ6": (1, 6), "ϵ16": (1, 6)}}
    """
    def __init__(self, s_n):
        super(EInfo, self).__init__()
        self.table_type = "ϵ_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": dict,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()

    def _init_cgc_static_db_folder(self):
        super(EInfo, self)._init_db_folder()
        # 初始化Sn=1时的ϵ_dict
        from core.cgc_utils.cgc_local_db import get_ϵ_file_name

        _, file_name = get_ϵ_file_name(1, [1], [1], [1], None)
        _, is_exist = self.exist_by_file_name(file_name)
        if is_exist is False:
            ϵ_dict = {"ϵ1": 1, "ϵ4": 1, "ϵ14": 1, "ϵ5": 1, "ϵ15": 1, "ϵ6": 1, "ϵ16": 1}
            table = {"file_name": file_name,
                     "data": ϵ_dict,
                     "flags": {}}
            flag, msg = self.insert(table)
            if not flag:
                raise Exception(msg)
            flag, msg = self.insert_txt(table)
            if not flag:
                raise Exception(msg)


class ISFInfo(CGCLocalDb):
    """
    这个db用来存ISF

    <CG>/isf_info/Sn/[σ]_[μ]/[ν’].pkl
    {
    "file_name": "Sn/[σ]_[μ]/[ν’]",
    "data": isf_square_dict,
    "flags": {"speed_time": speed_time}
    }

    isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                       "cols": [[ν], ([ν], β), ...],                   # 有自由度tuple，无自由度list
                       "isf": isf_square_matrix}                       # np.array([len(rows), len(cols)], dtype=float)

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；[ν]表示内积后的可能构型； (Sn)
    [σ’][μ’]表示参与内积的两个置换群构型；[ν’]表示内积后的可能构型； (Sn-1)
    beta和beta'分别对应[ν]以及[ν’]的多重性
    isf_square_dict数值是isf的平方，符号是平方前isf系数的符号
    len(rows) = len(cols)
    speed_time表示计算用时（秒）

    例如，
    <CG>/isf_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 1].pkl
    {"rows": [([3,1],[3,1]), ([3,1],[2,1,1]), ([2,1,1],[3,1]), ([2,1,1],[2,1,1])],
     "cols": [[4,1], ([3,2],1), ([3,2],2), [3,1,1]],
     "isf": sp.Matrix([[5/12, 1/2, 1/12, 0],
                       [-1/12, 0, 5/12, 1/2],
                       [-1/12, 0, 5/12, -1/2],
                       [5/12, -1/2, 1/12, 0]])}

    正则ISF索引的全部参数为：σ σ' μ μ' ν β ν' β'
    表示：|σ σ'> * |μ μ'> 的结果中，|νβ ν'β'>，的ISF系数平方
    """
    def __init__(self, s_n):
        super(ISFInfo, self).__init__()
        self.table_type = "isf_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": dict,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class EigenvaluesInfo(CGCLocalDb):
    """
    这个db用来存二循环类的本征值

    二循环类的本征值的落盘格式为：
    <CG>/eigenvalues_info/Sn.pkl         ->
    {
    "file_name": "Sn",
    "data": eigenvalues_list,             # list(int)
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    eigenvalues_list是yd按照Yamanouchi序排列对应的二循环类的本征值；
    speed_time表示计算用时（秒）

    例如：
    <CG>/eigenvalues_info/S6.pkl
    [15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15]

    另存：
    <CG>/eigenvalues_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": [],
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n},
              "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典
    """
    def __init__(self, s_n):
        super(EigenvaluesInfo, self).__init__()
        self.table_type = "eigenvalues_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": list,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class CGSeriesInfo(CGCLocalDb):
    """
    这个db用来存CG序列

    CG序列的落盘格式为：
    <CG>/cg_series_info/Sn/[σ]_[μ].pkl        ->
    {
    "file_name": "Sn/[σ]_[μ]",
    "data": cg_series,
    "flags": {"speed_time": speed_time}
    }

    其中，
    Sn表示n阶置换群;
    [σ][μ]表示参与内积的两个置换群构型；
    [ν]表示内积后的可能构型；
    cg_series是[ν]的线性组合系数；
    speed_time表示计算用时（秒）

    例如：
    <CG>/cg_series_info/Sn/[3]_[2, 1].pkl
    [0, 1, 0]

    另存：
    <CG>/cg_series_info/Finish_Sn.pkl ->
    {
    "file_name": "Finish_Sn",
    "data": np.array([0]),
    "flags": {"finish_s_n": s_n,
              "history_times": {"S{n}": time_of_s_n},
              "young_diagram_index": "young diagram list of Sn by young-yamanouchi"}
    }

    其中，
    Finish_Sn是固定名称，所有计算结果里都有
    s_n表示当前计算完成的最大阶数，是int型
    time_of_s_n表示Sn的计算时间，int型，外部是不断update的字典
    """

    def __init__(self, s_n):
        super(CGSeriesInfo, self).__init__()
        self.table_type = "cg_series_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": np.ndarray,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class CharacterAndGiInfo(CGCLocalDb):
    """
    这个db用来存特征标矩阵和与之对应的g_i

    它们的落盘格式为：
    <CG>/characters_and_gi_info/Sn.pkl       ->
    {
    "file_name": "Sn",
    "data": {"character": character_matrix, "gi": gi_array},
    "flags": {"speed_time": speed_time,
              "young_diagram_index": young_diagram_list_by_yamanouchi}  # 作为参考，非数据
    }

    其中，
    Sn表示n阶置换群;
    character_matrix是按照Littlewood书中给的表格和Yamanouchi序存放特征标矩阵;  # np.ndarray(int)
    gi是已经和特征标矩阵的列对齐的gi，同样来自Littlewood的书中;                 # np.ndarray(int)
    speed_time表示计算用时（秒）

    例如：
    <CG>/characters_and_gi_info/S4.pkl:
    {"character": [[ 1  1  1  1  1]
                   [ 3  1  0 -1 -1]
                   [ 2  0 -1  0  2]
                   [ 3 -1  0  1 -1]
                   [ 1 -1  1 -1  1]],
     "gi": [1, 6, 8, 6, 3]}

    另存：
    <CG>/characters_and_gi_info/Finish_Sn.pkl ->
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
        super(CharacterAndGiInfo, self).__init__()
        self.table_type = "characters_and_gi_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": dict,
            "flags": dict
        })
        self.s_n = s_n
        self.txt_limit = default_s_n
        self._init_cgc_static_db_folder()


class YamanouchiMatrixInfo(CGCLocalDb):
    """
    TODO tmp 转int型分数 考虑用fractions.Fraction

    这个db用来存Yamanouchi对换矩阵（仅限(ij)和(in)）

    Yamanouchi对换矩阵的落盘格式为：
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/ij(i,j).pkl 或
    <CG>/yamanouchi_matrix_info/Sn/[ν_i]/in(i,n).pkl     ->
    {
    "file_name": "Sn/[ν_i]/ij(i,j)",
    "data": matrix_ij 或 matrix_in,         #np.ndarray(float)    #它的维度可以在<CG>/young_tableaux_info/Sn/[ν_i]_num里取得
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
            "data": (list, np.ndarray, sp.Matrix),
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
    "data": total_num
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    total_num就是len({"m_i": young_table})，表示杨盘总数;

    例如：
    S3/[2, 1]_num.pkl: 2

    另存：杨盘相位因子
    <CG>/young_tableaux_info/Sn/[ν_i]_Λ.pkl  ->
    {
    "file_name": "Sn/[ν_i]_Λ",
    "data": phase_factor_list
    "flags": {}
    }

    其中，
    Sn表示n阶置换群;
    [ν_i]表示杨图;
    phase_factor_list就是Λ按照m从小到大，或者说Yamanouchi序排列的

    例如：
    S3/[2, 1]_Λ.pkl: [1, -1]
    """

    def __init__(self, s_n):
        super(YoungTableInfo, self).__init__()
        self.table_type = "young_tableaux_info"
        self.map_id = "file_name"
        self.design_table_type.update({
            "data": (dict, int, list),
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
