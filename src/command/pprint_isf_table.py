# -*- coding: utf-8 -*-
"""
打印可读性较好的ISF
p.s.
python pprint_isf_table.py --Sn 2 --σ 2 --μ 2 --ν_st 1
python pprint_isf_table.py --Sn 3 --σ 2 1 --μ 2 1 --ν_st 2
python pprint_isf_table.py --Sn 6 --σ 3 2 1 --μ 3 2 1 --ν_st 4 1
"""


import argparse
from pprint import pprint
import sympy as sp
from sympy import Rational as Ra
from sympy import sqrt
from core.isf_and_cgc import load_isf


def main(args):
    # 防bug
    '''
    这个bug来自sympy1.0版本之后的类注册器。如果不手动建立一次非整数的Rational，让它的注册表里存在Rational，直接读pickle会报错。

    C, including its class ClassRegistry, has been deprecated since SymPy
    1.0. It will be last supported in SymPy version 1.0. Use direct
    imports from the defining module instead. See
    https://github.com/sympy/sympy/issues/9371 for more info.
    '''
    _ = sp.Matrix([[Ra(1) / 2]])

    isf_param = (args.Sn, args.σ, args.μ, args.ν_st)  # TODO 参数检查
    flag, isf_square_dict = load_isf(*isf_param)
    if flag is False or not isinstance(isf_square_dict, dict):
        print("Load isf Fail: {}".format(isf_square_dict))
        return
    # 后面假设数据正确，就不一一验证了
    rows = isf_square_dict.get("rows")
    cols = isf_square_dict.get("cols")
    isf = isf_square_dict.get("isf")
    # sp.Matrix支持str输出，所以矩阵扩列就可以了
    str_rows = sp.Matrix([str(i) for i in rows])
    str_cols = sp.Matrix([["ISF"] + [str(i) for i in cols]])
    isf_insert_1_col = isf.col_insert(0, str_rows)  # 行的名字，加在列上
    isf_insert_1_col_1_row = isf_insert_1_col.row_insert(0, str_cols)  # 列的名字，加在行上
    pprint(isf_insert_1_col_1_row)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--Sn',
        type=int,
        help='Sn of ISF',
        required=True
    )
    parser.add_argument(
        '--σ',
        type=int,
        nargs="+",
        help='σ of ISF',
        required=True
    )
    parser.add_argument(
        '--μ',
        type=int,
        nargs="+",
        help='μ of ISF',
        required=True
    )
    parser.add_argument(
        '--ν_st',
        type=int,
        nargs="+",
        help='ν_st of ISF',
        required=True
    )
    args = parser.parse_args()
    main(args)
