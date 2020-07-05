import os
import logging
from plumbum import local
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from plumbum.machines.local import LocalCommand
from EukMetaSanity.src.utils.helpers import touch
from dask.distributed import Client, wait, as_completed
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""


class OutputResultsFileError(FileNotFoundError):
    pass


class Task(ABC):
    def __init__(self, input_path_dict: Dict[str, List[str]], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, mode: int, required_data: List[str]):
        # Store passed input flag:input_path dict
        # Require input file name passed
        for data in required_data:
            assert data in input_path_dict.keys(), data
        self._name = db_name
        self._input_path_dict = input_path_dict
        # Instantiate output dict variable
        self.required_data = [Data.OUT]
        self._output_paths_dict: Dict[str, List[str]] = {}
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, ConfigManager.THREADS))
        # Store path manager
        self._pm = pm
        # Store config manager
        self._cfg = cfg
        # Store primary calling program(s)
        self._prog = cfg.config.get(db_name, ConfigManager.PATH)
        self._prog2 = None
        self._prog3 = None
        self._prog4 = None
        for _path, _attr in (
            (ConfigManager.PATH2, "_prog2"),
            (ConfigManager.PATH3, "_prog3"),
            (ConfigManager.PATH4, "_prog4"),
        ):
            if _path in cfg.config[db_name].keys():
                setattr(self, _attr, cfg.config.get(db_name, _path))
        # Developer(0) or User(1) mode
        self._mode = mode
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        # Store id of record in Task
        self._record_id = record_id
        super().__init__()

    @property
    def name(self):
        return self._name

    @property
    def output(self):
        return self._output_paths_dict

    @output.setter
    def output(self, v: Dict[str, List[str]]):
        self._output_paths_dict = v

    @property
    def program(self):
        return local[self._prog]

    @property
    def program2(self):
        return local[self._prog2]

    @property
    def program3(self):
        return local[self._prog3]

    @property
    def program4(self):
        return local[self._prog4]

    @property
    def input(self) -> Dict[str, List[str]]:
        return self._input_path_dict

    @property
    def mode(self) -> int:
        return self._mode

    @property
    def record_id(self) -> str:
        return self._record_id

    @property
    def cfg(self) -> ConfigManager:
        return self._cfg

    @property
    def wdir(self) -> str:
        return self._wdir

    @property
    def threads(self) -> int:
        return self._threads_pw

    @property
    def pm(self) -> PathManager:
        return self._pm

    @abstractmethod
    def run(self) -> None:
        # Ensure that required_data is set
        assert self.required_data is not None
        for data in self.required_data:
            # Alert for missing required data output
            assert data in self._output_paths_dict.keys(), "Missing required %s" % data
        # Check if task has completed based on provided output data
        completed = True
        for _path in self._output_paths_dict[Data.OUT]:
            if not os.path.exists(_path):
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

    @abstractmethod
    def results(self) -> Dict[str, List[str]]:
        # Check that all required datasets are fulfilled
        for data in self.required_data:
            # Alert for missing required data output
            assert data in self._output_paths_dict.keys(), "Missing required %s" % data
            # Alert if data output is provided, but does not exist
            for _path in self._output_paths_dict[data]:
                if not os.path.exists(_path):
                    # Write dummy file if in developer mode
                    if self._mode == 0:
                        touch(_path)
                    else:
                        raise OutputResultsFileError(_path)
        return self._output_paths_dict

    @staticmethod
    # Function logs and runs dask command
    def log_and_run(cmd: LocalCommand, test: int):
        logging.info(str(cmd))
        if test == 1:
            cmd()


class TaskList(ABC):
    def __init__(self, task_list: List[Task], statement: str, workers: int, cfg: ConfigManager, pm: PathManager,
                 mode: int):
        self._statement = statement
        self._tasks: List[Task] = task_list
        self._workers = workers
        self._cfg = cfg
        self._pm = pm
        # Single(0) or Threaded(1) mode
        self._mode = mode
        super().__init__()

    @property
    def cfg(self):
        return self._cfg

    @property
    def pm(self):
        return self._pm

    @property
    def tasks(self):
        return self._tasks

    @abstractmethod
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
            for _task in as_completed(futures):
                _task.result()
            wait(futures)
            client.close()

    @abstractmethod
    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        # Run task list
        results = (task.results() for task in self._tasks)
        return (
            [result[Data.OUT] for result in results],  # Output files using required Data object
            self.cfg,
            self.pm,
            [task.record_id for task in self.tasks],
            self._mode,
        )


if __name__ == "__main__":
    pass
