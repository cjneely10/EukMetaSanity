import os
import logging
from plumbum import local
from abc import ABC, abstractmethod
from EukMetaSanity.src.utils.data import Data
from typing import Dict, List, Tuple, Callable
from plumbum.machines.local import LocalCommand
from EukMetaSanity.src.utils.helpers import touch
from dask.distributed import Client, wait, as_completed
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""


def program_catch(f: Callable):
    def _add_try_except(self, *args, **kwargs):
        try:
            f(self, *args, **kwargs)
        except ProcessExecutionError as e:
            logging.info(e)
            print(e)
    return _add_try_except


class OutputResultsFileError(FileNotFoundError):
    pass


class Task(ABC):
    def __init__(self, input_path_dict: Dict[Data.Type, List[str]], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, mode: int, required_data: List[str]):
        # Store passed input flag:input_path dict
        # Require input file name passed
        for data in required_data:
            assert data in input_path_dict.keys(), data
        self._name = db_name
        self._input_path_dict = input_path_dict
        # Instantiate output dict variable
        self._required_data = [Data.Type.OUT]
        self._output_paths_dict: Dict[Data.Type, List[str]] = {}
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
    def added_flags(self) -> List[str]:
        return self.cfg.get_added_flags(self.name)

    @property
    def name(self) -> str:
        return self._name

    @property
    def output(self) -> Dict[Data.Type, List[str]]:
        return self._output_paths_dict

    @output.setter
    def output(self, v: Dict[Data.Type, List[str]]):
        self._output_paths_dict = v

    @property
    def program(self) -> LocalCommand:
        return local[self._prog]

    @property
    def program2(self) -> LocalCommand:
        return local[self._prog2]

    @property
    def program3(self) -> LocalCommand:
        return local[self._prog3]

    @property
    def program4(self) -> LocalCommand:
        return local[self._prog4]

    @property
    def input(self) -> Dict[Data.Type, List[str]]:
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

    def results(self) -> Dict[Data.Type, List[str]]:
        # Check that all required datasets are fulfilled
        for data in self._required_data:
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

    # Function logs and runs dask command
    def log_and_run(self, cmd: LocalCommand):
        logging.info(str(cmd))
        if self.mode == 1:
            cmd()

    @abstractmethod
    def run(self) -> None:
        # Ensure that required_data is set
        assert self._required_data is not None
        for data in self._required_data:
            # Alert for missing required data output
            assert data in self._output_paths_dict.keys(), "Missing required %s" % data
        # Check if task has completed based on provided output data
        completed = True
        for _path in self._output_paths_dict[Data.Type.OUT]:
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
    def __init__(self, new_task: type, input_paths: List[List[str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str], mode: int, data_function: Callable,
                 required_data: Dict[Data.Type, List[str]] = None):
        if required_data is None:
            required_data = {}
        # Call data function for pertinent info
        name, _, statement = data_function()
        # Get workers for TaskList
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        # Get log statement
        self._statement = statement % (workers, int(cfg.config.get(name, ConfigManager.THREADS)))
        # Store list of tasks to complete
        self._tasks: List[Task] = [
            new_task(
                {Data.Type.IN: (input_path if isinstance(input_path, list) else [input_path]), **required_data},
                cfg,
                pm,
                record_id,
                name,
                mode,
                required_data=[Data.Type.IN, *required_data.keys()],
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

    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        # Run task list
        results = (task.results() for task in self._tasks)
        return (
            [result[Data.Type.OUT] for result in results],  # Output files using required Data object
            self.cfg,
            self.pm,
            [task.record_id for task in self.tasks],
            self._mode,
        )


if __name__ == "__main__":
    pass
