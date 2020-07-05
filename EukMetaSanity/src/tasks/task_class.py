import os
import logging
from plumbum import local
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from dask.distributed import Client, wait, as_completed
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""

touch = local["touch"]


class OutputResultsFileError(FileNotFoundError):
    pass


class Task(ABC):
    def __init__(self, input_path_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, required_data: List[str], mode: int):
        # Store passed input flag:input_path dict
        # Require input file name passed
        for data in required_data:
            assert data in input_path_dict.keys(), data
        self._input_path_dict = input_path_dict
        # Instantiate output dict variable
        self.required_data = [Data.OUT]
        self.output_paths_dict: Dict[str, str] = {}
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, "THREADS"))
        # Store path manager
        self._pm = pm
        # Store config manager
        self._cfg = cfg
        # Developer(0) or User(1) mode
        self._mode = mode
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        self._record_id = record_id
        super().__init__()

    @property
    def input(self):
        return self._input_path_dict

    @property
    def mode(self):
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
    def results(self) -> Dict[str, str]:
        # Check that all required datasets are fulfilled
        for data in self.required_data:
            # Alert for missing required data output
            assert data in self.output_paths_dict.keys(), "Missing required %s" % data
            # Alert if data output is provided, but does not exist
            # Write dummy file if in developer mode
            if not os.path.exists(self.output_paths_dict[data]):
                if self._mode == 0:
                    touch[self.output_paths_dict[data]]()
                else:
                    raise OutputResultsFileError(self.output_paths_dict[data])
        return self.output_paths_dict


class TaskList(ABC):
    def __init__(self, task_list: List[Task], statement: str, workers: int, cfg: ConfigManager, pm: PathManager,
                 mode: int):
        self._tasks: List[Task] = task_list
        logging.info(statement)
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
        # # Single
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
    def results(self):
        # Gather results to list and return
        return [task.results() for task in self._tasks]

    @abstractmethod
    def output(self) -> Tuple[List[str], ConfigManager, PathManager, List[str], int]:
        # Run task list
        return (
            [result[Data.OUT] for result in self.results()],  # Output files using required Data object
            self.cfg,
            self.pm,
            [task.record_id for task in self.tasks],
            self._mode,
        )


if __name__ == "__main__":
    pass
