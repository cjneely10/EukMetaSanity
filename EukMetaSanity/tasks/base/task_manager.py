import os
from shutil import copy
from typing import List, Dict
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.config_manager import ConfigManager
from EukMetaSanity.tasks.base.dependency_graph import DependencyGraph
from EukMetaSanity.tasks.manager.pipeline_manager import PipelineManager


class TaskManager:
    def __init__(self, pm: PipelineManager, cfg: ConfigManager, pam: PathManager,
                 input_files: List[Dict[str, Dict[str, object]]], input_prefixes: List[str], debug: bool, command: str):
        self.dep_graph = DependencyGraph(pm.programs[command])
        self.task_list = self.dep_graph.sorted_tasks()
        self.completed_tasks: List[TaskList] = []
        self.pm = pam
        self.cfg = cfg
        self.debug = debug
        self.command = command
        self.input_files = input_files
        self.input_prefixes = input_prefixes

    def run(self, output_dir: str):
        task = self.task_list[0](self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
        task.run()
        self.completed_tasks.append(task)
        i = 1
        while i < len(self.task_list):
            task = self.task_list[i](self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
            for req_str in task.requires:
                task.update({req_str: self.dep_graph.idx[req_str].output})
            task.run()
            self.completed_tasks.append(task)
        self.summarize(os.path.join(output_dir, "results", self.command), self.command)

    def summarize(self, _final_output_dir: str, _name: str):
        if not os.path.exists(_final_output_dir):
            os.makedirs(_final_output_dir)
        _paths_output_file = open(os.path.join(os.path.dirname(_final_output_dir), "%s-paths_summary.tsv" % _name), "w")
        for task_list in self.completed_tasks:
            for task_result, task_record_id in zip(task_list.output()):
                if "final" in task_result.keys():
                    _sub_out = os.path.join(_final_output_dir, task_record_id)
                    if not os.path.exists(_sub_out):
                        os.makedirs(_sub_out)
                    for _file in task_result["final"]:
                        # Write info to file
                        if isinstance(_file, dict) and len(_file.keys()) > 0:
                            _sorted_keys = sorted(list(_file.keys()))
                            # Header
                            _paths_output_file.write("".join(("\t".join(["ID"] + _sorted_keys), "\n")))
                            # Path info
                            _paths_output_file.write(
                                "".join((
                                    "\t".join((
                                        task_record_id,  # Name of record
                                        *(os.path.join(_sub_out, str(_file[_f])) for _f in _sorted_keys)
                                    )), "\n"
                                ))
                            )
                        # Generate link of path, or create copy as requested
                        elif isinstance(_file, str):
                            if os.path.exists(_file):
                                copy(_file, _sub_out)
