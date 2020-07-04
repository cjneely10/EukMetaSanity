from typing import Dict, List
from abc import ABC, abstractmethod
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager


class Task(ABC):
    def __init__(self, input_path_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str):
        # Store passed input flag:input_path dict
        self._input_path_dict = input_path_dict
        # Instantiate output dict variable
        self.output_paths_dict: Dict[str, str] = {}
        # Store threads and workers
        self._threads_pw = int(cfg.config.get(db_name, "THREADS")),
        self._workers = int(cfg.config.get(db_name, "WORKERS")),
        # Store path manager
        self._pm = pm
        # Store config manager
        self._cfg = cfg
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        super().__init__()

    @property
    def cfg(self):
        return self._cfg

    @property
    def wdir(self):
        return self._wdir

    @property
    def threads(self):
        return self._threads_pw

    @property
    def workers(self):
        return self._workers

    @property
    def pm(self):
        return self._pm

    @abstractmethod
    def run(self):
        assert self.output_paths_dict != {}

    @abstractmethod
    def results(self):
        if self.output_paths_dict == {}:
            self.run()
        return self.output_paths_dict


class TaskList(ABC):
    def __init__(self, task_list: List[Task]):
        self._tasks: List[Task] = task_list
        super().__init__()

    @property
    def tasks(self):
        return self._tasks

    @abstractmethod
    def run(self):
        for task in self._tasks:
            task.run()


if __name__ == "__main__":
    pass
