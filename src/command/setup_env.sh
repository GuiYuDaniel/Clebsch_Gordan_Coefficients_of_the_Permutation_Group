#! /bin/bash
# -*- coding: utf-8 -*-

#=====================================================
# this shell for:
#  1, setup conda env
#  2, setup PYTHONPATH
#  3, cd to the command directory
# to copy this shell to your home folder or any other easy seen folder is helpful
#=====================================================


### params
_PROJECT_DIR="/Users/guiyu/PycharmProjects"  # change to your project dir !
_PROJECT_NAME="Clebsch_Gordan_Coefficients_of_the_Permutation_Group"  # if the project is renamed, change to the new name
_WORKING_DIR_NAME="test"  # an existed dir in project, p.s.[test, src/command, src, ...], the shell will cd to this dir
_CONDA_ENV_NAME="cgc"  # change to your conda env name !


### shell
echo "shell for:"
echo ""

echo "1, setup ${_CONDA_ENV_NAME} envs!"
conda activate ${_CONDA_ENV_NAME}
echo ""

echo "2, cd to working dir:"
cd ${_PROJECT_DIR}/${_PROJECT_NAME}/${_WORKING_DIR_NAME}
echo `pwd`
echo ""

echo "3, set PYTHONPATH as:"
_PYTHONPATH_1="${_PROJECT_DIR}/${_PROJECT_NAME}"  # top path
_PYTHONPATH_2="${_PROJECT_DIR}/${_PROJECT_NAME}/src"  # src path
_PYTHONPATH_3="${_PROJECT_DIR}/${_PROJECT_NAME}/test"  # test path
export PYTHONPATH=${_PYTHONPATH_1}:${_PYTHONPATH_2}:${_PYTHONPATH_3}
echo $PYTHONPATH
echo ""
