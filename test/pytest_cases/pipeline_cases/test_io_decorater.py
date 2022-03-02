# -*- coding: utf-8 -*-


from functools import wraps


def _check_and_change_params_another(func, param_2="default_param_2", param_3="default_param_3"):

    param_list = [param_2, param_3]
    if any(not isinstance(param_i, str) for param_i in param_list):
        param_type_list = [str(type(i)) for i in param_list]
        err_msg = "param_list={} must all str but not with param_type_list={}".format(param_list, param_type_list)
        return False, err_msg
    if not callable(func):
        err_msg = "function name={} is not callable, pls check".format(func.__name__)
        return False, err_msg

    @wraps(func)
    def wrapper(param_1):
        if any(not isinstance(param_i, str) for param_i in [param_1, param_2, param_3]):
            err_msg = "{} {} {} not all str".format(param_1, param_2, param_3)
            return False, err_msg
        param_1 = param_1 + "_and_" + param_2 + "_and_" + param_3
        result = func(param_1, param_2, param_3)
        return result

    return wrapper


class DecoraterExample(object):
    """
    这个装饰器包括：
    1，载入默认参数
    2，提前返回
    3，修改参数
    """

    # 这是个装饰器，按照python指南，应该写在类的外部，但是从内容上考量，写在里面更易读，所以这里稍微违背一下python指南
    # TODO 写装饰器staticmethod会导致不合法，暂时只想到这种直接不加的写法
    def _check_and_change_params(func):
        @wraps(func)
        def wrapper(param_1, param_2=func.__defaults__[0], param_3=func.__defaults__[1]):
            if any(not isinstance(param_i, str) for param_i in [param_1, param_2, param_3]):
                err_msg = "{} {} {} not all str".format(param_1, param_2, param_3)
                return False, err_msg
            param_1 = param_1 + "_and_" + param_2 + "_and_" + param_3
            result = func(param_1, param_2, param_3)
            return result
        return wrapper

    @staticmethod
    @_check_and_change_params
    def print_params(param_1, param_2="default_param_2", param_3="default_param_3"):
        return param_1, param_2


class TestDecoraterExample(object):

    def setup_class(self):

        self.user_param_1 = "user_param_1"
        self.user_param_2 = "user_param_2"
        self.user_param_3 = "user_param_3"
        self.default_param_2 = "default_param_2"
        self.default_param_3 = "default_param_3"
        self.answer_1 = "user_param_1_and_default_param_2_and_default_param_3"
        self.answer_12 = "user_param_1_and_user_param_2_and_default_param_3"
        self.answer_13 = "user_param_1_and_default_param_2_and_user_param_3"
        self.answer_123 = "user_param_1_and_user_param_2_and_user_param_3"
        self.wrong_param = None
        self.wrong_answer = False

    def teardown_class(self):
        pass

    def test_with_user_params_1(self):
        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1)
        assert rst_1 == self.answer_1
        assert rst_2 == self.default_param_2

    def test_with_user_params_12(self):
        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1, self.user_param_2)
        assert rst_1 == self.answer_12
        assert rst_2 == self.user_param_2

        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1, param_2=self.user_param_2)
        assert rst_1 == self.answer_12
        assert rst_2 == self.user_param_2

    def test_with_user_params_13(self):
        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1, param_3=self.user_param_3)
        assert rst_1 == self.answer_13
        assert rst_2 == self.default_param_2

    def test_with_user_params_123(self):
        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1, self.user_param_2, self.user_param_3)
        assert rst_1 == self.answer_123
        assert rst_2 == self.user_param_2

        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1, self.user_param_2, param_3=self.user_param_3)
        assert rst_1 == self.answer_123
        assert rst_2 == self.user_param_2

        rst_1, rst_2 = DecoraterExample.print_params(self.user_param_1,
                                                     param_2=self.user_param_2,
                                                     param_3=self.user_param_3)
        assert rst_1 == self.answer_123
        assert rst_2 == self.user_param_2

    def test_with_wrong_user_params(self):
        rst_1, rst_2 = DecoraterExample.print_params(self.wrong_param)
        assert rst_1 == self.wrong_answer
        assert isinstance(rst_2, str)


class TestIODecorater(object):
    """test decoraters in io.py"""

    def setup_class(self):
        pass

    def teardown_class(self):
        pass
