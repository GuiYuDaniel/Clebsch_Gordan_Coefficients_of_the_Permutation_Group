# -*- coding: utf-8 -*-


import time


def get_run_function_n_times(func, *args, n=100, **kwargs):
    """获取run n次function的时间，内部函数，为了时间更加准确，暂不检查参数了"""
    start_time = time.time()
    for i in range(int(n)):
        func_rst = func(*args, **kwargs)
    speed_time = time.time() - start_time
    return "func={} runs {} times with speeding {}s".format(func.__name__, n, speed_time), func_rst
