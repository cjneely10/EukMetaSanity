import os
import json
from shutil import copy
from typing import List, Dict
from collections import defaultdict
from EukMetaSanity.tasks.utils.helpers import touch
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.config_manager import ConfigManager
from EukMetaSanity.tasks.base.dependency_graph import DependencyGraph
from EukMetaSanity.tasks.manager.pipeline_manager import PipelineManager


class TaskManager:
    """ Class interfaces stored pipelines and input data from users.

    Dependency graph created using pipeline TaskList class object `requires` members, and this graph is
    output in topologically-sorted order to run pipeline

    Output of each task is automatically parsed into other tasks that depend on its output

    At the end of a pipeline, any task that contains a "final" key in its `output` member will have valid member
    paths copied into a final results directory.

    All dictionary data will also be serialized into a final "task.json" file for loading into other pipelines

    """
    def __init__(self, pm: PipelineManager, cfg: ConfigManager, pam: PathManager,
                 input_files: List[Dict[str, Dict[str, object]]], input_prefixes: List[str], debug: bool, command: str):
        self.dep_graph = DependencyGraph(pm.programs[command])
        self.task_list = self.dep_graph.sorted_tasks()
        self.completed_tasks: Dict[str, TaskList] = {}
        self.pm = pam
        self.cfg = cfg
        self.debug = debug
        self.command = command
        self.input_files = input_files
        self.input_prefixes = input_prefixes

    def run(self, output_dir: str):
        """ Run pipeline!

        :param output_dir: Directory to write final result file
        """
        task = self.task_list[0](self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
        task.run()
        self.completed_tasks[task.name] = task
        i = 1
        while i < len(self.task_list):
            task = self.task_list[i](self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
            for req_str in task.requires:
                for out_data in self.completed_tasks[req_str].output()[0]:
                    task.update({req_str: out_data})
            task.run()
            self.completed_tasks[task.name] = task
            i += 1
        self.summarize(os.path.join(output_dir, "results", self.command), self.command)

    def summarize(self, _final_output_dir: str, _name: str):
        """ Summarize contents of pipeline - populate output files and summary .json paths file

        :param _final_output_dir: Output directory
        :param _name: Name of pipeline
        """
        if not os.path.exists(_final_output_dir):
            os.makedirs(_final_output_dir)
        _paths_output_file = open(os.path.join(os.path.dirname(_final_output_dir), "%s.json" % _name), "w")
        output_dicts: Dict[str, Dict[str, object]] = defaultdict(default_factory=defaultdict(str))
        # Collect header ids and corresponding files
        for task_list in self.completed_tasks.values():
            output = task_list.output()
            for task_result, task_record_id in zip(*output):
                if "final" in task_result.keys():
                    _sub_out = os.path.join(_final_output_dir, task_record_id)
                    if not os.path.exists(_sub_out):
                        os.makedirs(_sub_out)
                    for _file in task_result["final"]:
                        output_dicts[task_record_id] = {task_list.name: {_file: task_result[_file]}}
                        _file = task_result[_file]
                        if isinstance(_file, str):
                            if os.path.exists(_file):
                                copy(_file, _sub_out)
                            elif self.debug:
                                touch(os.path.join(_sub_out, _file))
        del output_dicts["default_factory"]
        # Also store original files
        for in_prefix, in_data in zip(self.input_prefixes, self.input_files):
            output_dicts[in_prefix].update(in_data["root"])
        json.dump(dict(output_dicts), _paths_output_file, sort_keys=True)
        _paths_output_file.close()