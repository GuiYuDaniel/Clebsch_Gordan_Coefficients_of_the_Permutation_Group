说明：

计算结果：<top>/cgc_results  简记为：<CG>

计算程序：<top>/src

工具：python36

结果说明：(详见：<top>/src/core/cgc_utils/cgc_db_typing.py)(结果的计算顺序从下至上)

    CG系数存放于：C:\CG\...
        无β：../CG/CG/Sn/[σ]_[μ]/[ν]m.pkl         # {(m1,m2):C,’N’:N}

    ISF存放于：../CG/ISF/Sn/[σ]_[μ]／[ν]—->[ν’]/[σ’]_[μ’].pkl

        ../CG/ISF/Sn/[σ]_[μ]／[ν]_[ν’]/file_list.pkl     # 记录全部有效[[σ’],[μ’]]

        ../CG/ISF/Sn/[σ]_[μ]／folder_list.pkl              # 记录全部有效[[ν],[ν’]]


    CGorder存放于：../CG/CG_order/Sn/[σ]_[μ]/[ν].pkl      # /Users/guiyu/Desktop/CG/CG_order/S3/[2,1]_[2,1]/[3].pkl

    特征标表存放于：../CG/KAI/Sn/[r].pkl
               ../CG/KAI/Sn/[r].txt(手动录入的，便于检查)
               ../CG/KAI/Sn/gi.pkl
               ../CG/KAI/Sn/gi.txt（手动录入的，便于检查）
        (P205,《群论习题精解》,马中骐；按列存入,编程检查正交归一)

    二循环类算符本征值存放于：/Users/guiyu/Desktop/CG/eigenvalue/Sn/[ri].pkl    #如/Users/guiyu/Desktop/CG/eigenvalue/S3/[2，1].pkl
                                                                                      #存放S3群[2,1]构型的二循环类算符本征值
        同目录下：/Users/guiyu/Desktop/CG/eigenvalue/Sn/[ri].txt                        #同内容txt文件

    本征值to构型存放于：/Users/guiyu/Desktop/CG/eigenvalue2r/Sn/λ.pkl    #如/Users/guiyu/Desktop/CG/eigenvalue2r/S3/0.pkl
                                                                                      #存放S3群二循环类算符本征值等于0的构型————[2,1]
        同目录下：/Users/guiyu/Desktop/CG/eigenvalue2r/Sn/λ.txt                         #同内容txt文件

    对换矩阵存放于：
        (i,j)对换：C:\CG\YamanouchiMatrix(ij)\Sn\ij[i,j]\[ri].pkl    #如C:\CG\YamanouchiMatrix(ij)\S4\ij[2,3]\[2,1,1].pkl
                                                                                       #存放S4群[2,1,1]构型的（23）交换矩阵
            同目录下C:\CG\YamanouchiMatrix(ij)\Sn\ij[i,j]\[ri].npy                      #是在python中易于存取的相同文件
            同目录下C:\CG\YamanouchiMatrix(ij)\Sn\ij[i,j]\[ri].txt                      #存放Sn<=5的Sn群[ri]构型的（i，j）交换矩阵

        (i,n)对换：C:\CG\YamanouchiMatrix(in)\Sn\in[i,n]\[ri].pkl    #如C:\CG\YamanouchiMatrix(in)\S4\in[2,4]\[2,1,1].pkl
                                                                                      #存放S4群[2,1,1]构型的（24）交换矩阵
            同目录下C:\CG\YamanouchiMatrix(in)\Sn\in[i,n]\[ri].npy                      #是在python中易于存取的相同文件
            同目录下C:\CG\YamanouchiMatrix(in)\Sn\in[i,n]\[ri].txt                      #存放Sn<=5的Sn群[ri]构型的（i，n）交换矩阵
                                                                                      #存放S4群[2,1,1]构型的（24）交换矩阵

    对换矩阵存放于：
        (i,j)对换：<CG>/yamanouchi_matrix_info/Sn/[ν_i]/ij(i,j).pkl #如<CG>/yamanouchi_matrix_info/S4/[2,1,1]/ij(2,3).pkl
            数据结构：{"data": matrix_ij}           #list(list(float))                  #存放S4群[2,1,1]构型的（23）交换矩阵
            txt展示：value of "data"
            #[[-0.5, 0.8660254037844386, 0], [0.8660254037844386, 0.5, 0], [0, 0, -1]]

        (i,n)对换：<CG>/yamanouchi_matrix_info/Sn/[ν_i]/in(i,n).pkl #如<CG>/yamanouchi_matrix_info/S4/[2,1,1]/in(2,4).pkl
            数据结构：{"data": matrix_in}           #list(list(float))                  #存放S4群[2,1,1]构型的（24）交换矩阵
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

