说明：

计算结果：<top>/cgc_results  简记为：<CG>

计算程序：<top>/src

工具：python36

结果说明：(详见：<top>/src/core/cgc_utils/cgc_db_typing.py)(结果的计算顺序从下至上)
# TODO numpy转sympy的例子都要修改

    # TODO 矩阵优化存储
    # "从S7的[4,2,1]_[4,2,1] [3,2,1,1]_τ8_m2看 35*35=1225个组合中击中510个float
    #  字典需要18528B 而35*35维np仅仅需要9912(无论填什么) "
    CG系数存放于：（对应表4.13）
        <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_τ_m.pkl  # 无多重性，则为[ν]_m.pkl
        # {(σ_m, μ_m): cgc_square_value}
        # 其中，cgc_square_value是cg系数的平方乘cg系数的符号，并且，它们的绝对值之和应该等于1
        数据结构：{"data": cgc_square_dict}            #dict(tuple(int, int): fraction)
        txt展示：value of "data"
        #<CG>/cgc_info/S4/[2, 2]_[3, 1]/[2, 1, 1]_m2.pkl
        #{(2, 1): 1/2, (1, 3): 1/4, (2, 2): 1/4}
        #<CG>/cgc_info/S5/[3, 1, 1]_[3, 2]/[1, 1, 1]_τ1_m4.pkl
        #{(1, 4): 1/8, (2, 4): -1/16, (3, 5): 1/16, (4, 1): -1/10, (4, 2): -1/5,
          (4, 4): 3/80, (5, 3): 1/5, (5, 5): -3/80, (6, 3): -1/10, (6, 5): -3/40}

    ISF存放于：（对应表4.19）
        <CG>/isf_info/Sn/[σ]_[μ]/[ν’].pkl
        # <CG>/isf_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 1].pkl
        数据结构{"data": {"rows": row_list,
                        "cols": col_list,
                        "isf": isf_square_dict}}
        isf_square_dict = {"rows": [([σ'], [μ'], τ'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                           "cols": [[ν], ([ν], τ), ...],                   # 有自由度元组，无自由度列表, len(rows)=len(cols)
                           "isf": isf_square_matrix}                       # sp.Matrix([len(rows)]) fraction
        # {"data": {"rows": [([3,1],[3,1]), ([3,1],[2,1,1]), ([2,1,1],[3,1]), ([2,1,1],[2,1,1])],
                    "cols": [[4,1], ([3,2],1), ([3,2],2), [3,1,1]],
                    "isf": sp.Matrix([[5/12, 1/2, 1/12, 0],
                                      [-1/12, 0, 5/12, 1/2],
                                      [-1/12, 0, 5/12, -1/2],
                                      [5/12, -1/2, 1/12, 0]])}}

    ϵ存放于：
        <CG>/ϵ_info/Sn/[σ]_[μ]/[ν]_τ.pkl  # 无多重性，则为[ν].pkl
        # <CG>/ϵ_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 2]_2.pkl
        数据结构{"data": ϵ_dict}
        ϵ_dict = {"ϵx": ϵx: int}  # 主数据，记录ϵx的值
        flags = {"ϵx": (m_σ, m_μ,): tuple}  # 辅助数据，记录对应CGC的来源，便于计算对称性ϵ_dict
        # [3, 1, 1]_[3, 1, 1]/[3, 2]_2
        {"data": {"ϵ0": 1, "ϵ1": 1, "ϵ4": -1, "ϵ14": -1,
                  "ϵ5": 1, "ϵ15": -1, "ϵ6": -1, "ϵ16": 1},
         "flags": {"ϵ0": (1, 1), "ϵ1": (1, 1), "ϵ4": (6, 6), "ϵ14": (6, 6),
                   "ϵ5": (6, 4), "ϵ15": (3, 1), "ϵ6": (1, 3), "ϵ16": (4, 6)}}
        # [3, 1, 1]_[3, 1, 1]/[2, 1, 1]_2
        {"data": {"ϵ0": 1, "ϵ1": -1, "ϵ4": -1, "ϵ14": 1,
                  "ϵ5": -1, "ϵ15": -1, "ϵ6": 1, "ϵ16": 1},
         "flags": {"ϵ0": (1, 4), "ϵ1": (4, 1), "ϵ4": (6, 3), "ϵ14": (3, 6),
                   "ϵ5": (6, 1), "ϵ15": (6, 1), "ϵ6": (1, 6), "ϵ16": (1, 6)}}

    元组合与对称组合存放于：
        <CG>/symmetry_info/Sn/meta_σμν.pkl
        数据结构{"data": meta_σμν_tuple_list}
        meta_σμν_tuple_list = [(meta_σ, meta_μ, meta_ν), ...]
        例如：
        # <CG>/symmetry_info/S3/meta_σμν.pkl
        # {"data": [([3], [3], [3]), ([3], [2, 1], [2, 1]), ([2, 1], [2, 1], [2, 1])]}

        <CG>/symmetry_info/Sn/[σ][μ][ν]_symmetries.pkl
        数据结构{"data": sym_σμν_dict}
        sym_σμν_dict = {(j_σ_s, j_μ_s, j_ν_s): [meta_key]}  # 来源可以不唯一
        例如：
        # <CG>/symmetry_info/S3/[2, 1][2, 1][2, 1]_symmetries.pkl
        # {"data": {('[2, 1]', '[2, 1]', '[2, 1]'): ["σμν", "μσν", ...]}}  # 这个例子里，24个meta_key都是

    二循环类算符本征值存放于：(本征值to构型不存了，使用函数计算)
        <CG>/eigenvalues_info/Sn.pkl                #<CG>/eigenvalues_info/S6.pkl
        数据结构：{"data": eigenvalues}               #list(int)
        txt展示：value of "data"
        #[15, 9, 5, 3, 3, 0, -3, -3, -5, -9, -15]

    CG序列存放于：
        <CG>/cg_series_info/Sn/[σ]_[μ].pkl           #<CG>/cg_series_info/Sn/[3]_[2, 1].pkl
        数据结构：{"data": cg_series_array}            #np.ndarray(int)
        txt展示：value of "data"
        #[0, 1, 0]

    特征标表和gi存放于：
        <CG>/characters_and_gi_info/Sn.pkl          #如<CG>/characters_and_gi_info/S4.pkl
        #按照Littlewood书中给的表格和Yamanouchi序存放特征标矩阵和gi
        数据结构：{"data":
                  {"character": character_matrix,   #np.ndarray(int)    #特征标矩阵
                   "gi": gi_array}}                 #np.ndarray(int)    #已经和特征标矩阵的列对齐的gi
        # flags中会列出young_diagram_index
        txt展示：value of "data"
        #{"character": [[ 1  1  1  1  1]
                        [ 3  1  0 -1 -1]
                        [ 2  0 -1  0  2]
                        [ 3 -1  0  1 -1]
                        [ 1 -1  1 -1  1]],
          "gi": [1, 6, 8, 6, 3]}

    对换矩阵存放于：
        (i,j)对换：<CG>/yamanouchi_matrix_info/Sn/[ν_i]/ij(i,j).pkl #如<CG>/yamanouchi_matrix_info/S4/[2,1,1]/ij(2,3).pkl
            数据结构：{"data": matrix_ij}           #np.ndarray(float)                  #存放S4群[2,1,1]构型的（23）交换矩阵
            txt展示：value of "data"
            #[[-0.5, 0.8660254037844386, 0], [0.8660254037844386, 0.5, 0], [0, 0, -1]]

        (i,n)对换：<CG>/yamanouchi_matrix_info/Sn/[ν_i]/in(i,n).pkl #如<CG>/yamanouchi_matrix_info/S4/[2,1,1]/in(2,4).pkl
            数据结构：{"data": matrix_in}           #np.ndarray(float)                  #存放S4群[2,1,1]构型的（24）交换矩阵
            txt展示：value of "data"
            #[[-0.5, 0.28867513459481287, -0.8164965809277259],
              [0.28867513459481287, -0.8333333333333334, -0.4714045207910317],
              [-0.8164965809277259, -0.4714045207910317, 0.3333333333333333]]
        注：有且仅有 ij(Sn-1,Sn) 与 in(Sn-1,Sn) 是重复的

    杨盘存放于：<CG>/young_tableaux_info/Sn/[ν_i].pkl  #如<CG>/young_tableaux_info/S3/[2,1].pkl  #存放S3群[2,1]构型的所有杨盘
        数据结构：{"data": {"m_i": young_table}}       #"2":[[1,3],[2]](str:list(list(int)))    #杨盘编号：杨盘
        txt展示：value of "data"
        #{"1": [[1, 2], [3]], "2": [[1, 3], [2]]}

        另存：杨盘总数/也是不可约表示维度h/也是最大m编号
        <CG>/young_tableaux_info/Sn/[ν_i]_num.pkl #如<CG>/young_tableaux_info/S3/[2,1]_num.pkl #存放S3群[2,1]构型的杨盘总数
        数据结构：{"data": total_num}                  #2(int)
        txt展示：value of "data"
        #2
        #len({"m_i": young_table}) = total_num

        另存：杨盘相位因子
        <CG>/young_tableaux_info/Sn/[ν_i]_Λ.pkl  #如<CG>/young_tableaux_info/S3/[2,1]_Λ.pkl  #存放S3群[2,1]构型的相位因子
        数据结构：{"data": phase_factor_list}           #如[1, -1]
        txt展示：value of "data"

    分支律存放于：<CG>/branching_laws_info/Sn/[ν_i].pkl   #如<CG>/branching_laws_info/S4/[2,2].pkl  #存放S4群[2,2]构型的分支律
        数据结构：{"data": {"BL_num": bl_num,         #1(int)                     #S4群[2,2]构型分支律的分支数
                          "rows": rows,             #[1](list(int))             #S4群[2,2]构型分支律去掉格子的py行号列表
                          "cols": cols,             #[1](list(int))             #S4群[2,2]构型分支律去掉格子的py列号列表
                          "before_YD": [ν’]}}       #[[2, 1]](list(list(int)))  #S4群[2,2]构型分支律前置杨盘列表
        txt展示：{"BL_num": bl_num, "rows": rows, "cols": cols, "before_YD": [ν’]}

    杨图存放于：<CG>/young_diagrams_info/S{n}.pkl          #如<CG>/young_diagrams_info/S4.pkl           #存放S4群的所有构型
        数据结构：{"data": [ν]}                           #[[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]](list(list(int)))
        txt展示：[ν]


思路说明：
    ISF方法由 Cn-1 造 Cn

    Core:1，由[σ][μ]根据CG序列定[ν],τ可能，再由[σ][μ][ν],τ定ISF计算目标
         2，由[ν]根据分支律定[ν’]，list，再由[σ][μ]根据分支律定[σ‘][μ’]list
         3，由(in),λn,λn-1，求ISF:[ν],τ,[ν’],τ’,[σ][σ’],[μ][μ’]
               处理(λn)-(λn-1)撞衫
         4，由[ν][σ][μ]定m，m1，m2
         5，由ISF，CGn-1:[ν’]m’,τ’,[σ’]m1’,[μ’]m2’，求 CGn:[ν]m,τ,[σ]m1,[μ]m2

    亮点：1，相位的继承分支律

程序说明：

