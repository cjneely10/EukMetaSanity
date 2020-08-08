import os
import logging
from plumbum import local, BG
from abc import ABC, abstractmethod
from dask.distributed import Client, wait
from typing import Dict, List, Tuple, Callable
from EukMetaSanity.tasks.manager.data import Data
from EukMetaSanity.tasks.utils.helpers import touch
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.config_manager import ConfigManager
from plumbum.commands.processes import ProcessExecutionError
from plumbum.machines.local import LocalCommand, LocalMachine

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""


def program_catch(f: Callable):
    def _add_try_except(self, *args, **kwargs):
        try:
            return f(self, *args, **kwargs)
        except ProcessExecutionError as e:
            logging.info(e)
            print(e)
    return _add_try_except


class OutputResultsFileError(FileNotFoundError):
    pass


class Task(ABC):
    def __init__(self, input_path_list: List[str], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, mode: int):
        # Store passed input flag:input_path dict
        self._name = db_name
        self._input_path_list = input_path_list
        # Instantiate output dict variable
        self._output_paths: List[object] = []
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, ConfigManager.THREADS))
        # Store path manager
        self._pm = pm
        # Store config manager
        self._cfg = cfg
        # Dynamically generate program attributes for easy in API access
        self._set_api_accessors(cfg, db_name)
        # Developer(0) or User(1) mode
        self._mode = mode
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        # Store id of record in Task
        self._record_id = record_id
        self.skip = False
        super().__init__()

    def _set_api_accessors(self, cfg: ConfigManager, db_name: str):
        for _path, _value in cfg.config[db_name].items():
            # Set attribute
            if _path.split("_")[0].isupper() and _value != "None":
                _path = _path.lower()
                # Set as self.program is only default Config PATH variable
                # Set attribute for ease of use in API
                _set_attr = _value
                # Add program from local environment
                if _path.startswith("program"):
                    # assert os.path.exists(_value)
                    _set_attr = local[_set_attr]
                # Check for existence if a data value
                elif _path.startswith("data"):
                    for _val in _value.split(","):
                        assert os.path.exists(_val)
                setattr(
                    self,
                    _path,  # Name: PATH -> program/data; PATH2/DATA_2 = program2/data_2;
                    _set_attr,  # Local path, or config path, for calling program
                )

    @property
    def local(self) -> LocalMachine:
        return local

    @property
    def added_flags(self) -> List[str]:
        return self.cfg.get_added_flags(self.name)

    @property
    def name(self) -> str:
        return self._name

    @property
    def output(self) -> List[object]:
        return self._output_paths

    @output.setter
    def output(self, v: List[object]):
        self._output_paths = v

    @property
    def input(self) -> List[str]:
        return self._input_path_list

    @property
    def mode(self) -> int:
        return self._mode

    @property
    def record_id(self) -> str:
        return self._record_id

    @property
    def cfg(self) -> ConfigManager:
        return self._cfg

    # Returns dict of config section for this task
    @property
    def config(self) -> Dict[str, str]:
        return self._cfg.config[self._name]

    @property
    def wdir(self) -> str:
        return self._wdir

    @property
    def pm(self) -> PathManager:
        return self._pm

    def results(self) -> List[object]:
        # Check that all required datasets are fulfilled
        # Alert if data output is provided, but does not exist
        for _path in self._output_paths:
            if self.skip == "False" and isinstance(_path, str) and not os.path.exists(_path):
                # Write dummy file if in developer mode
                if self._mode == 0:
                    touch(_path)
                else:
                    raise OutputResultsFileError(_path)
        return self._output_paths

    # Function logs and runs dask command
    def log_and_run(self, cmd: LocalCommand):
        print("  " + str(cmd))
        logging.info(str(cmd))
        if self.mode == 1:
            logging.info(cmd())

    def batch(self, cmds: List[LocalCommand]):
        for i in range(0, len(cmds), int(self.threads)):
            running = []
            for j in range(i, i + int(self.threads)):
                if j >= len(cmds):
                    break
                print("  " + str(cmds[j]))
                logging.info(str(cmds[j]))
                if self.mode == 1:
                    f = cmds[j] & BG
                    running.append(f)
            all([_f.wait() for _f in running])

    @abstractmethod
    def run(self) -> None:
        # Check if task has completed based on provided output data
        completed = True
        for _path in self._output_paths:
            if isinstance(_path, str) and not os.path.exists(_path):
                # Only call function if missing path
                # Then move on
                completed = False
                break
        # Run if not completed (e.g. missing data)
        if completed:
            logging.info("%s  %s is complete" % (self.record_id, self.name))
        else:
            logging.info("%s  Running %s" % (self.record_id, self.name))
            # Gather all functions of the form run_1, run_2, etc.
            runnables = sorted(
                [func for func in dir(self) if func.startswith("run_")],
                key=lambda _f: int(_f.split("_")[1]),
            )
            # Run each function in series
            for func in runnables:
                if func.startswith("run_"):
                    getattr(self, func)()


class TaskList(ABC):
    def __init__(self, new_task: type, name: str, cfg: ConfigManager, input_paths: List[List[str]], pm: PathManager,
                 record_ids: List[str], mode: int):
        # Call data function for pertinent info
        dt = Data(cfg, name)
        name, _, statement = getattr(dt, dt.name)()
        # Get workers for TaskList
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        # Get log statement
        self._statement = statement % (workers, int(cfg.config.get(name, ConfigManager.THREADS)))
        # Store list of tasks to complete
        self._tasks: List[Task] = [
            new_task(
                (input_path if isinstance(input_path, list) else [input_path]),
                cfg,
                pm,
                record_id,
                name,
                mode,
            )
            for input_path, record_id in zip(input_paths, record_ids)
        ]
        # Store workers
        self._workers = workers
        # Store ConfigManager object
        self._cfg = cfg
        # Store PathManager object
        self._pm = pm
        # Single(0) or Threaded(1) mode
        self._mode = mode

    @property
    def cfg(self):
        return self._cfg

    @property
    def pm(self):
        return self._pm

    @property
    def tasks(self):
        return self._tasks

    def run(self):
        # Single
        logging.info(self._statement)
        print(self._statement)
        if self._mode == 0:
            for task in self._tasks:
                task.run()
        # Threaded
        else:
            futures = []
            client = Client(n_workers=self._workers, threads_per_worker=1)
            # Run each future
            for _task in self._tasks:
                futures.append(client.submit(_task.run))
            wait(futures)
            client.close()

    def output(self) -> Tuple[ConfigManager, List[List[object]], PathManager, List[str], int]:
        # Run task list
        return (
            self.cfg,
            [task.results() for task in self._tasks],  # Output files using required Data object
            self.pm,
            [task.record_id for task in self._tasks],
            self._mode,
        )

    # Write summary file of results
    def summarize(self, _final_output_dir: str, _name: str):
        # Create softlinks (or copies) of final output files to output directory
        _output = self.output()
        _output_files_list = _output[1]
        _files_prefixes = _output[3]
        if not os.path.exists(_final_output_dir):
            os.makedirs(_final_output_dir)
        _paths_output_file = open(os.path.join(os.path.dirname(_final_output_dir), "%s-paths_summary.tsv" % _name), "w")
        for _files, _file_prefix, _task in zip(_output_files_list, _files_prefixes, self.tasks):
            # Create subdirectory
            _sub_out = os.path.join(_final_output_dir, _file_prefix)
            if not os.path.exists(_sub_out):
                os.makedirs(_sub_out)
            # Copy results to results dir for easier access
            for _file in _files:
                # Write info to file
                if isinstance(_file, dict):
                    _sorted_keys = sorted(list(_file.keys()))
                    # Header
                    _paths_output_file.write("".join(("\t".join(["ID"] + _sorted_keys), "\n")))
                    # Path info
                    _paths_output_file.write(
                        "".join((
                            "\t".join((
                                _file_prefix,  # Name of record
                                *(os.path.join(_sub_out, str(_file[_f])) for _f in _sorted_keys)  # Files for record
                            )), "\n"
                        ))
                    )
                # Generate link of path, or create copy as requested
                elif isinstance(_file, str):
                    if os.path.exists(_file):
                        _task.log_and_run(
                            _task.program[
                                (*_task.added_flags),
                                _file, _sub_out,
                            ]
                        )
        _paths_output_file.close()


if __name__ == "__main__":
    pass
