说明：

计算结果：<top>/cgc_results  简记为：<CG>

计算程序：<top>/src

工具：python36

结果说明：(详见：<top>/src/core/cgc_utils/cgc_db_typing.py)(结果的计算顺序从下至上)

    # TODO tmp!
    # TODO [σ]_[μ] 与 [μ]_[σ]
    # TODO 尝试整分数化（可能不行，会遭遇多重根）
    # TODO 因为行列是关于σ_m, μ_m的全排列，而且数值化的时候，都有值，也就是说是不是用矩阵表示更好？
    # "从S7的[4,2,1]_[4,2,1] [3,2,1,1]_8_m2看 35*35=1225个组合中击中510个float
    #  字典需要18528B 而35*35维np仅仅需要9912(无论填什么) "
    CG系数存放于：
        <CG>/cgc_info/Sn/[σ]_[μ]/[ν]_β_m.pkl      # {(σ_m, μ_m): cgc_value, N: 平方和}   # TODO 矩阵优化存储
        #<CG>/cgc_info/S4/[2,2]_[3,1]/[2,1,1]_1_m2.pkl
        数据结构：{"data": cgc_square_dict}            #dict(tuple(int, int): float; str: float)
        txt展示：value of "data"
        #{(2, 1): 0.4999999999999999, (1, 3): 0.24999999999999994, (2, 2): 0.24999999999999994, 'N': 0.9999999999999998}

    # TODO tmp!
    # TODO [σ]_[μ] 与 [μ]_[σ]
    # TODO 尝试整分数化
    ISF存放于：
        <CG>/isf_info/Sn/[σ]_[μ]/[ν’].pkl  # 存形如表4.19那样的矩阵
        # <CG>/isf_info/S5/[3, 1, 1]_[3, 1, 1]/[3, 1].pkl  # 表4.19-17b
        数据结构{"data": {"rows": row_list,
                        "cols": col_list,
                        "isf": isf_square_dict}}
        isf_square_dict = {"rows": [([σ'], [μ'], β'), ([σ'], [μ']), ...],  # 有自由度len3，无自由度len2
                           "cols": [[ν›], ([ν], β), ...],                   # 有自由度tuple，无自由度list
                           "isf": isf_square_matrix}                       # np.array([len(rows), len(cols)], dtype=float)
        # {"data": {"rows": [([3,1],[3,1]), ([3,1],[2,1,1]), ([2,1,1],[3,1]), ([2,1,1],[2,1,1])],
                    "cols": [[4,1], ([3,2],1), ([3,2],2), [3,1,1]],
                    "isf": np.array([[5/12, 1/2, 1/12, 0],
                                     [-1/12, 0, 5/12, 1/2],
                                     [-1/12, 0, 5/12, -1/2],
                                     [5/12, -1/2, 1/12, 0]])}}
        # <CG>/isf_info/S6

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

        另存：杨盘总数
        <CG>/young_tableaux_info/Sn/[ν_i]_num.pkl #如<CG>/young_tableaux_info/S3/[2,1]_num.pkl #存放S3群[2,1]构型的杨盘总数
        数据结构：{"data": {"total_num": total_num}}                  #2(int)
        txt展示：value of "data"
        #{"total_num": 2}

        #len({"m_i": young_table}) = total_num

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

    Core:1，由[σ][μ]根据CG序列定[ν],β可能，再由[σ][μ][ν],β定ISF计算目标
         2，由[ν]根据分支律定[ν’]，list，再由[σ][μ]根据分支律定[σ‘][μ’]list
         3，由(in),λn,λn-1，求ISF:[ν],β,[ν’],β’,[σ][σ’],[μ][μ’]
               处理(λn)-(λn-1)撞衫
         4，由[ν][σ][μ]定m，m1，m2
         5，由ISF，CGn-1:[ν’]m’,β’,[σ’]m1’,[μ’]m2’，求 CGn:[ν]m,β,[σ]m1,[μ]m2

    亮点：1，相位的继承分支律

程序说明：


φ

