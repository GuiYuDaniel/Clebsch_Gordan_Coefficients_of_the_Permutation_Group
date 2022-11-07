# -*- coding:utf8 -*-
"""
register some cgc config here
using these config as global params by python import
"""


import os
from itertools import product
from utils.config import singleton_config


"""
calculate params

这些参数是和计算CG系数有关的默认参数，
推荐的方式是在命令行中根据计算需要，输入非默认参数
用户也可以在这里修改默认参数，修改前请确保知晓修改参数的数学意义
"""
# Sn 置换群阶数
default_s_n = 7

# 符号计算/数值计算标志
default_symbolic_calculation = True


"""
system params

无需修改这些参数，系统就可以正常执行
如果希望修改，请确保用户知晓修改参数的意义，对当前程序、历史结果的影响
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
young_tableaux_phase_factor_file_name_format = os.path.join("S{}", "{}_Λ")  # "S{}/{}_Λ"
# yamanouchi matrix file name format
yamanouchi_matrix_file_name_format = os.path.join("S{}", "{}", "{}{}")  # "S{}/{}/{}{}"
# characters and gi file name format
characters_and_gi_file_name_format = "S{}"
# CG series file name format
cg_series_file_name_format = os.path.join("S{}", "{}_{}")  # "S{}/{}_{}"
# eigenvalues file name format
eigenvalues_file_name_format = "S{}"
# isf file name format
isf_file_name_format = os.path.join("S{}", "{}_{}", "{}\'")
# ϵ file name format
ϵ_file_name_format_none = os.path.join("S{}", "{}_{}", "{}")
ϵ_file_name_format = os.path.join("S{}", "{}_{}", "{}_τ{}")
# cgc file name format
cgc_file_name_format_none = os.path.join("S{}", "{}_{}", "{}_m{}")
cgc_file_name_format = os.path.join("S{}", "{}_{}", "{}_τ{}_m{}")
# symmetry combination file name format
meta_σμν_file_name_format = os.path.join("S{}", "meta_σμν")
sym_σμν_file_name_format = os.path.join("S{}", "{}{}{}_symmetries")

# min Sn 各计算单元的最小合法Sn
min_s_n_of_young_diagram = 1
min_s_n_of_branching_law = 1  # 理论上是2; 添加1是一种"解析延拓"
min_s_n_of_young_table = 1
min_s_n_of_yamanouchi_matrix = 2
min_s_n_of_characters_and_gi = 1
min_s_n_of_cg_series = 1
min_s_n_of_eigenvalue = 1  # 理论上是2; 添加1是一种"解析延拓"
min_s_n_of_isf = 2
min_s_n_of_ϵ = 2  # 尽管有1，但它是预置的，不是计算得来的
min_s_n_of_cgc = 1
min_s_n_of_sym = 1

# 2*2*2*6/2=24个ϵ的四分量映射
'''
其中， 
行将按照D3对称群命名[e, d, f , a, b, c](abc也对应了书中的123)
列保持书中的命名方式[0, 4, 5, 6]
以b4为例，就是先按照b的方式将σμν调整为νμσ，再按4的方式将νμσ前两个杨图取共轭得到ν~μ~σ
四元tuple为其可操作表示，
tuple[0]为σμν次序，后面的三个01数字表示按照当前顺序下，每个杨盘是自身0，还是共轭1
特别标记ϵ的是为了与书中记法对应

(e,0,0,0)    (a,0,0,0)    (b,0,0,0)    (c,0,0,0)    (d,0,0,0)    (f,0,0,0)  
 σμν: e0      μσν: a0      νμσ: b0      σνμ: c0      νσμ: d0      μνσ: f0
              ϵ1           ϵ2           ϵ3
 
(e,1,1,0)    (a,1,1,0)    (b,1,1,0)    (c,1,1,0)    (d,1,1,0)    (f,1,1,0)  
 σ~μ~ν:e4     μ~σ~ν:a4     ν~μ~σ:b4     σ~ν~μ:c4     ν~σ~μ:d4     μ~ν~σ:f4
 ϵ4
 
(e,1,0,1)    (a,1,0,1)    (b,1,0,1)    (c,1,0,1)    (d,1,0,1)    (f,1,0,1)  
 σ~μν~:e5     μ~σν~:a5     ν~μσ~:b5     σ~νμ~:c5     ν~σμ~:d5     μ~νσ~:f5
 ϵ5
 
(e,0,1,1)    (a,0,1,1)    (b,0,1,1)    (c,0,1,1)    (d,0,1,1)    (f,0,1,1)  
 σμ~ν~:e6     μσ~ν~:a6     νμ~σ~:b6     σν~μ~:c6     νσ~μ~:d6     μν~σ~:f6
 ϵ6
'''
# group_d3 = ("σμν", "μσν", "νμσ", "σνμ", "νσμ", "μνσ")
# group_d3 * group_k4就可以造24种ϵ了
group_d3 = ((0, 1, 2), (1, 0, 2), (2, 1, 0), (0, 2, 1), (2, 0, 1), (1, 2, 0))  # 0表示σ；1表示μ；2表示ν
group_k4 = ((False, False, False), (True, True, False), (True, False, True), (False, True, True))
