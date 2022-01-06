# -*- coding: utf-8 -*-


import os
import json
from utils.log import get_logger


logger = get_logger(__name__)


class Json(object):
    """
    file to json, this json could with comments which start with #
    """

    @staticmethod
    def file_to_json(file_path):
        """
        Read a file and try to transfer to json format
        :param file_path: file full path(string)
        :return: json data from file(dict)
        """
        if not isinstance(file_path, str):
            err_msg = "file_path={} must be string".format(file_path)
            logger.error(err_msg)
            raise Exception(err_msg)
        if os.path.exists(file_path) is False or os.path.isfile(file_path) is False:
            err_msg = "file_path={} not exist or not file, pls check".format(file_path)
            logger.error(err_msg)
            raise Exception(err_msg)
        with open(file_path, 'rt', encoding='UTF-8') as f:
            data = f.read()
            try:
                json_data = json.loads(data)
            except json.decoder.JSONDecodeError as e:
                logger.exception(e)
                raise e
        return json_data

    @classmethod
    def _remove_comment(cls, json_data):
        """
        Remove 2 kinds of situation:
        1, {"### remove ###": ..., "keep": ...} ==> {"keep": ...}
        2, ["### remove ###", "keep", keep]     ==> ["keep", keep]

        Yah,
        hopeful comments should be added only in LIST and KEY of dict,
        not in string or value of dict.

        :param json_data: original json file data(dict list str)
        :return: pure json data(dict list str)
        """
        if isinstance(json_data, dict):
            all_keys = list(json_data.keys())
            for key in all_keys:
                if isinstance(key, str) and key.startswith("#"):
                    # delete1: key-value which key has comment
                    del json_data[key]
                    continue

                val = json_data.get(key)
                if isinstance(val, dict):
                    json_data[key] = cls._remove_comment(val)
                elif isinstance(val, list):
                    json_data[key] = cls._remove_comment(val)
            return json_data

        elif isinstance(json_data, list):
            rst_list = []
            for item in json_data:
                if isinstance(item, str):
                    if item.startswith("#"):
                        continue
                    else:
                        rst_list.append(item)

                if isinstance(item, dict):
                    rst_list.append(cls._remove_comment(item))
                elif isinstance(item, list):
                    rst_list.append(cls._remove_comment(item))
            return rst_list

        else:
            return json_data

    @classmethod
    def file_to_json_without_comments(cls, file_path):
        json_data = cls.file_to_json(file_path)
        return cls._remove_comment(json_data)
