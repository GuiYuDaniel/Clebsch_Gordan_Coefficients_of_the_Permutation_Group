# -*- coding:utf8 -*-
"""
register some cgc config here
using these config as global params by python import
"""


import os
from utils.config import singleton_config


"""
calculate params

这些参数是和计算CG系数有关的默认参数，
推荐的方式是在命令行中根据计算需要，输入非默认参数
用户也可以在这里修改默认参数，修改前请确保知晓修改参数的数学意义
"""
# Sn 置换群阶数
default_s_n = 7


"""
system params

无需修改这些参数，系统就可以正常执行
如果希望修改，请确保用户知晓修改参数的意义
"""
# <CGC>
top_path = singleton_config.top_path
# <CGC>/cgc_results = <CGC_rst>
cgc_rst_folder = "cgc_results"
calculated_tables_txt_limit = -1  # -1表示永远存txt

# db table name
# young diagrams file name format
young_diagrams_file_name_format = "S{}"
# branching laws file name format
branching_laws_file_name_format = os.path.join("S{}", "{}")  # "S{}/{}"
# young tableaux file name format
young_tableaux_file_name_format = os.path.join("S{}", "{}")  # "S{}/{}"
young_tableaux_num_file_name_format = os.path.join("S{}", "{}_num")  # "S{}/{}_num"
# yamanouchi matrix file name format
yamanouchi_matrix_file_name_format = os.path.join("S{}", "{}", "{}{}")  # "S{}/{}/{}{}"
# characters and gi file name format
characters_and_gi_file_name_format = "S{}"
# CG series file name format
cg_series_file_name_format = os.path.join("S{}", "{}_{}")  # "S{}/{}_{}"
# eigenvalues file name format
eigenvalues_file_name_format = "S{}"
