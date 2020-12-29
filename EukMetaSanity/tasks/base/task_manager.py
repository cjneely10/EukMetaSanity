"""
Module contains TaskManager class
"""

import os
import pickle
from shutil import copy
from typing import List, Dict, Tuple, Iterable, Union
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.config_manager import ConfigManager
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
        self.task_list = self.dep_graph.sorted_tasks
        self.completed_tasks: Dict[Tuple[str, str], TaskList] = {}
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
        # Generate first task from class object
        task = self.task_list[0][0](
            self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug, self.task_list[0][1],
            [{} for _ in range(len(self.input_files))], [self.input_files[k]["root"] for k in range(len(self.input_files))])
        # Run and store results
        task.run()
        self.completed_tasks[(task.name, task.scope)] = task
        i = 1
        # Run each task in list
        while i < len(self.task_list):
            # Create reference to prior task
            old_task = task
            # Collect input data that requirements and dependencies request
            to_add = []
            # Collect `dependency_input` dict from completed task list
            expected_input = []
            for k in range(len(old_task.output()[1])):
                inner_add = {}
                # Requirements are stored at outermost scope
                for req_str in self.task_list[i][0].requires:
                    inner_add[req_str] = self.completed_tasks[(req_str, "")].tasks[k].output
                # Dependencies are stored at task scope level, but may also be from task's scope's scope
                for req_str in self.task_list[i][0].depends:
                    inner_add[req_str.name] = self.completed_tasks[
                        (req_str.name, task.scope) if (req_str.name, task.scope) in self.completed_tasks.keys()
                        else (req_str.name, task.name)
                    ].tasks[k].output
                to_add.append(inner_add)
                # Dependency input will either come from root or will be collected from a task that has already run
                if self.task_list[i][2] != "root":
                    expected_input.append(
                        self.completed_tasks[(self.task_list[i][2], "")].tasks[k].output
                    )
                else:
                    expected_input.append(self.input_files[k]["root"])
            # Generate next task based on input fromrequired dependencies/requirements
            task = self.task_list[i][0](
                self.cfg, self.input_files, self.pm, self.input_prefixes, self.debug, self.task_list[i][1],
                to_add, expected_input)
            # Run task and store object in completed task list
            task.run()
            self.completed_tasks[(task.name, task.scope)] = task
            i += 1
        # Summarize final results based on requested `final` labels
        self._summarize(os.path.join(output_dir, "results", self.command), self.command)

    def _summarize(self, _final_output_dir: str, _name: str):
        """ Summarize contents of pipeline - populate output files and summary .pkl file

        :param _final_output_dir: Output directory
        :param _name: Name of pipeline
        """
        # Create results directory
        if not os.path.exists(_final_output_dir):
            os.makedirs(_final_output_dir)
        # Collect items marked final from each completed task
        output_data: Dict[str, object] = {}
        for task_list in self.completed_tasks.values():
            output = task_list.output()
            i = 0
            for task_result, task_record_id in zip(*output):
                if "final" in task_result.keys():
                    # Copy output data and update output data dict to write
                    output_data.update(
                        self._manage_output(_final_output_dir, task_record_id, task_result, task_list.name, i)
                    )
                i += 1
        # Serialize output to file
        pickle.dump(output_data, open(os.path.join(_final_output_dir, self.command + ".pkl"), "wb"))

    def _manage_output(self, output_directory: str, record_id: str, task_result: Dict[str, Union[object, Iterable]],
                       task_name: str, completed_tasklist_idx: int) -> Dict[str, object]:
        """ Create output subdirectory and copy output marked `final` to output directory

        :param output_directory: Pipeline's output directory path
        :param record_id: ID of input data
        :param task_result: Output of task (from calling output())
        :param task_name: Name/class name assigned to task
        :param completed_tasklist_idx: Position in self.completed_tasks containing this record id's data
        :return:
        """
        final_output_paths: Dict[str, object] = {}
        _sub_out = os.path.join(output_directory, record_id)
        # Create task subdirectory in results directory
        if not os.path.exists(_sub_out):
            os.makedirs(_sub_out)
        for _file in task_result["final"]:
            # Check if file is output from current task
            if _file in task_result.keys():
                out_file = task_result[_file]
            # Otherwise file should exist from existing tasklist
            else:
                class_path = _file.split(".")
                try:
                    out_file = self.completed_tasks[
                        (".".join(class_path[0:-1]), task_name)
                    ].output()[0][completed_tasklist_idx][class_path[-1]]
                except KeyError as e:
                    print("Unable to locate task output %s" % _file)
                    # Raise fatal error if unable to handle final output request
                    raise RuntimeError from e
            # Store output file object
            final_output_paths[_file] = out_file
            # str object are checked to see if they are files and are written
            if isinstance(out_file, str):
                if os.path.exists(out_file):
                    copy(out_file, _sub_out)
        return final_output_paths

