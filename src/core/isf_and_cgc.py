# -*- coding: utf-8 -*-
"""
this code for creating ISF and CGC
"""

# 见《群表示论的新途径》陈金全（上海科学技术出版社1984）本段代码主要依据置换群CGC的第二种递推计算方法
# ISF：第四章第19节
# CGC：第四章第11、12节
# 整体思路是，先求Sn的ISF，再利用Sn的ISF以及Sn-1的CG，拼出Sn的CG


import copy
import time
from conf.cgc_config import default_s_n
from core.young_diagrams import load_young_diagrams, is_young_diagram
from core.cgc_utils.cgc_db_typing import EigenvaluesInfo
from core.cgc_utils.cgc_local_db import get_eigenvalues_file_name, get_eigenvalues_finish_s_n_name
from utils.log import get_logger


logger = get_logger(__name__)
