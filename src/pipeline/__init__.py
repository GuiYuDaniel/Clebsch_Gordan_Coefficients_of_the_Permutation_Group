# -*- coding:utf8 -*-
"""
pipeline负责作业的调度模块，
实现将DAG图指定的workflow，按照一定顺序（拓扑序）派发作业并执行的功能。

本模块目前是串行同步调度，视未来情况，可改为并行异步调度（celery，redis）

pipeline负责构建整体计算图
pipenode负责构建单个计算节点
pipetask负责执行以及调度

pipe构建中可以抛错并中断，但执行任务过程中，函数务必执行到底直到状态改变为成功、失败、中断等结束状态

调用时：
ppl = PipeLine(workflow_conf=conf, ppl_name=xxx, ppl_id=yyy)  # 创建ppl
ppt = PipeTask(ppl)  # 创建ppk
status, msg = ppt.start(inputs_of_first_node)  # 执行

重启：
ppt = PipeTask().load_by_id(ppt_id)  # 根据pipetask_id就可以恢复所有
status, msg = ppt.restart()  # 可重启的pipetask必须通过preparation阶段，此时一切运行必要信息都已计算和保存
"""
