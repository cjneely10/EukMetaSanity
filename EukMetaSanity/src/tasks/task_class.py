import logging
from typing import Dict, List
from abc import ABC, abstractmethod
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""


class Task(ABC):
    def __init__(self, input_path_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, required_data: List[str]):
        # Store passed input flag:input_path dict
        # Require input file name passed
        for data in required_data:
            assert data in input_path_dict.keys(), data
        self._input_path_dict = input_path_dict
        # Instantiate output dict variable
        self.required_data: List[str] = []
        self.output_paths_dict: Dict[str, str] = {}
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, "THREADS"))
        self._workers = int(cfg.config.get(db_name, "WORKERS"))
        # Store path manager
        self._pm = pm
        self._is_complete = False
        # Store config manager
        self._cfg = cfg
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        self._record_id = record_id
        super().__init__()

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
    def workers(self) -> int:
        return self._workers

    @property
    def pm(self) -> PathManager:
        return self._pm

    @abstractmethod
    def run(self) -> None:
        # Set complete once run is done
        self._is_complete = True

    @abstractmethod
    def results(self) -> Dict[str, str]:
        # Run if not done so already
        assert self._is_complete
        for data in self.required_data:
            assert data in self.output_paths_dict.keys(), data
        return self.output_paths_dict

    @abstractmethod
    def parse_output(self, output_files: List[str]) -> List[Dict[str, str]]:
        pass


class TaskList(ABC):
    def __init__(self, task_list: List[Task], statement: str):
        self._tasks: List[Task] = task_list
        logging.info(statement)
        super().__init__()

    @property
    def tasks(self):
        return self._tasks

    @abstractmethod
    def run(self):
        for task in self._tasks:
            task.run()

    @abstractmethod
    def results(self):
        # Gather results to list and return
        return [task.results() for task in self._tasks]


if __name__ == "__main__":
    pass
