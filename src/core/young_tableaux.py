# -*- coding: utf-8 -*-
"""
this code for creating Young Tableaux by Sn and young_diagrams
计算Sn杨图的杨盘
"""

# \subsection{杨盘}
#
# 在$S_n$群杨图$[\nu]$上按下面规则填上1，2，...，n 个数字所构成的一个图称为杨盘。
# 规则要求我们在任一行自左向右和任一列自上向下阅读时，这些数字必须是按增加次序排列的。
# 我们用$Y^{[\nu]}_m$来代表一个杨盘，m为杨盘的编号。
#
# 我们可用杨盘来标志$S_n$的标准基$|Y^{[\nu]}_m>$。
# 它的意义是按照这样来规定的：
# 在杨盘$Y^{[\nu]}_m$中去掉方块n，就得到涉及排列好的n-1个方块的一个盘$Y^{[\nu']}_{m'}$；
# 再去掉方块n-1，就又得到排列好的下一个盘$Y^{[\nu'']}_{m''}$，等等，直到最后。
# 规定不可约的基失$|Y^{[\nu]}_m>$是这样一种基失：
# 它分别属于群$S_n,S_{n-1},S_{n-2},...,S_1$的IR$[\nu],[\nu'],[\nu''],...,[1]$。


import time
from core.cgc_utils.cgc_db_typing import YoungDiagramInfo
from core.cgc_utils.cgc_local_db import get_young_diagrams_file_name, get_young_diagrams_finish_s_n_name
from conf.cgc_config import default_s_n
from utils.log import get_logger


logger = get_logger(__name__)


def calc_young_tableaux(s_n=default_s_n):
    pass