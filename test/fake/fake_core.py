# -*- coding:utf8 -*-
"""
模拟的core函数，用于测试pipeline，形式上这里的函数都应该有输入和输出(函数式函数，但可以改状态)
虽然pipeline调度不是靠输入输出确定顺序的，但，目前暂不设计非函数了（没想好支持更多的好处）
如果原始函数形式不是函数式，套个壳子嘛
"""
import os
import pwd
import time


def merge(param_1, param_2, flag=True):
    if flag:
        sleep_time = 1
        print("will sleep {}s".format(sleep_time))
        time.sleep(sleep_time)
        print("sleep {}s end".format(sleep_time))
        rst = str(param_1)+str(param_2)
    else:
        rst = False
    return rst


def self_defined_print(print_str, flag=True):
    if flag:
        sleep_time = 1
        print("will sleep {}s".format(sleep_time))
        time.sleep(sleep_time)
        print("sleep {}s end".format(sleep_time))
        print(print_str)
        return flag
    else:
        return False


def what_is_my_uid(flag=True):
    if flag:
        uid = os.getuid()
    else:
        uid = False
    return uid


def who_am_i(uid):
    if uid:
        try:
            pw_name = pwd.getpwuid(uid).pw_name
        except Exception as e:
            print(e)
            print("Cannot use uid to get user name, we will use another method")
            pw_name = pwd.getpwuid(os.getuid()).pw_name
    else:
        pw_name = pwd.getpwuid(os.getuid()).pw_name
    return pw_name


def where_am_i(normal_flag=True, point_path=None):
    """
    :param normal_flag:
    :param point_path:
    :return: the path of this function or point_path
    """
    if normal_flag:
        sleep_time = 1
        print("will sleep {}s".format(sleep_time))
        time.sleep(sleep_time)
        print("sleep {}s end".format(sleep_time))
        if not point_path:
            path = os.getcwd()
        else:
            path = os.path.abspath(point_path)
    else:
        raise Exception('Fake Error')
    return path


def what_is_the_time_now(flag=True):
    if flag:
        asc_time = time.asctime()
        sleep_time = 1
        print("will sleep {}s".format(sleep_time))
        time.sleep(sleep_time)
        print("sleep {}s end".format(sleep_time))
    else:
        return False
    return asc_time


def output_str_func():
    # 这个函数python里是可以给他传空tuple的
    # 这个函数就用无参数传入测，后面的测试还是用标准函数化的方式测
    return "this is str"


def is_input_str_func(string):
    return True if isinstance(string, str) else False


def output_all_str_func(flag):
    return "this is str 1", "this is str 2"


def is_input_all_str_func(*args):
    return True if all(isinstance(i, str) for i in args) else False
