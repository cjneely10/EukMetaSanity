import os
import logging
from plumbum import local
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Callable
from EukMetaSanity.src.utils.helpers import touch
from EukMetaSanity.src.tasks.manager.data import Data
from dask.distributed import Client, wait, as_completed
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from plumbum.machines.local import LocalCommand, LocalMachine
from EukMetaSanity.src.utils.config_manager import ConfigManager

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
        self._output_paths: List[str] = []
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
        super().__init__()

    def _set_api_accessors(self, cfg: ConfigManager, db_name: str):
        for _path, _value in cfg.config[db_name].items():
            # Set attribute
            if _path.split("_")[0].isupper() and _value != "None":
                _path = _path.lower()
                # Set as self.program is only default Config PATH variable
                # Set attribute for ease of use in API
                _set_attr = _value
                if _path.startswith("program"):
                    _set_attr = local[_set_attr]
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
    def output(self) -> List[str]:
        return self._output_paths

    @output.setter
    def output(self, v: List[str]):
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

    def results(self) -> List[str]:
        # Check that all required datasets are fulfilled
        # Alert if data output is provided, but does not exist
        for _path in self._output_paths:
            if not os.path.exists(_path):
                # Write dummy file if in developer mode
                if self._mode == 0:
                    touch(_path)
                else:
                    raise OutputResultsFileError(_path)
        return self._output_paths

    # Function logs and runs dask command
    def log_and_run(self, cmd: LocalCommand):
        logging.info(str(cmd))
        if self.mode == 1:
            cmd()

    @abstractmethod
    def run(self) -> None:
        # Check if task has completed based on provided output data
        completed = True
        for _path in self._output_paths:
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
            for _task in as_completed(futures):
                _task.result()
            wait(futures)
            client.close()

    def output(self) -> Tuple[ConfigManager, List[List[str]], PathManager, List[str], int]:
        # Run task list
        return (
            self.cfg,
            [task.results() for task in self._tasks],  # Output files using required Data object
            self.pm,
            [task.record_id for task in self.tasks],
            self._mode,
        )


if __name__ == "__main__":
    pass
