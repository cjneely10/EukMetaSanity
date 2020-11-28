import os
from typing import List, Dict
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.config_manager import ConfigManager
from EukMetaSanity.tasks.base.dependency_graph import DependencyGraph
from EukMetaSanity.tasks.manager.pipeline_manager import PipelineManager


class TaskManager:
    def __init__(self, pm: PipelineManager, cfg: ConfigManager, pam: PathManager,
                 input_files: List[Dict[str, Dict[str, object]]], input_prefixes: List[str], debug: bool, command: str):
        self.dep_graph = DependencyGraph(pm.programs[command])
        self.task_list = self.dep_graph.sorted_tasks().sort(reverse=True)
        self.pm = pam
        self.cfg = cfg
        self.debug = debug
        self.command = command
        self.input_files = input_files
        self.input_prefixes = input_prefixes

    def run(self, output_dir: str):
        task = self.task_list.pop()(self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
        task.run()
        while len(self.task_list) > 0:
            task = self.task_list.pop()(self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug)
            for req_str in task.requires:
                task.input.update({req_str: self.dep_graph.idx[req_str].output})
            task.run()
        task.summarize(os.path.join(output_dir, "results", self.command), self.command)
        # task = self.task_list[0](self.cfg, self.input_files[0], self.pm, self.input_prefixes[1], self.debug)
        # for i in range(1, len(self.task_list)):
        #     # Run task
        #     task.run()
        #     task = self.task_list[i](*task.output())
        # # Must call output on last task to generate final summary statistics
        # task.run()
        # # Create summary using final Summarize task
        # task.summarize(os.path.join(output_dir, "results", self.command), self.command)
