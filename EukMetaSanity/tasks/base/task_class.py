import os
import logging
from plumbum import local, BG
from abc import ABC, abstractmethod
from dask.distributed import Client, wait
from EukMetaSanity.tasks.helpers import touch
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.slurm_caller import SLURMCaller
from plumbum.machines.local import LocalCommand, LocalMachine
from typing import Dict, List, Tuple, Callable, Optional, Union, Iterable
from EukMetaSanity.tasks.base.config_manager import ConfigManager, MissingDataError

"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""


def program_catch(f: Callable):
    """ Decorator function checks for ProcessExecutionErrors and FileExistErrors when running an executable
    Logging info recorded and printed to stdout

    :param f: Function to test and whose exceptions to catch
    :return: Decorated function
    """

    def _add_try_except(self, *args, **kwargs):
        try:
            return f(self, *args, **kwargs)
        except ProcessExecutionError as e:
            logging.info(e)
            with open(os.path.join(self.wdir, "task.err"), "a") as w:
                w.write(str(e))
            print(e)
        except FileExistsError as e:
            logging.info(e)
            with open(os.path.join(self.wdir, "task.err"), "a") as w:
                w.write(str(e))
            print(e)
        except ValueError as e:
            logging.info(e)
            with open(os.path.join(self.wdir, "task.err"), "a") as w:
                w.write(str(e))
            print(e)

    return _add_try_except


class OutputResultsFileError(FileNotFoundError):
    pass


class Task(ABC):
    def __init__(self, input_data: Dict[str, Dict[str, object]], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, mode: int, scope: str):
        self._name = db_name
        # Store passed input flag:input_path dict
        self._input_data = input_data
        # Instantiate output dict variable
        self._output_paths: Dict[str, object] = {}
        # Store config manager
        self._cfg = cfg
        # Store threads and workers
        self._threads_pw = "1"
        self._scope = scope
        if db_name in self.cfg.config.keys():
            self._threads_pw = str(cfg.config[db_name][ConfigManager.THREADS])
        else:
            self._threads_pw = str(cfg.config[scope][ConfigManager.THREADS])
        # Store path manager
        self._pm = pm
        # Developer(0) or User(1) mode
        self._mode = mode
        # Add name of db
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        # Store id of record in Task
        self._record_id = record_id
        # Check if task is set to be skipped in config file
        self.is_skip = "skip" in self.config.keys() and self.config["skip"] is True
        super().__init__()

    @abstractmethod
    def run(self):
        pass

    def _run(self) -> None:
        # Check if task has completed based on provided output data
        if self.is_skip:
            return
        completed = True
        for _path in self._output_paths.values():
            if isinstance(_path, str) and not os.path.exists(_path):
                # Only call function if missing path
                # Then move on
                completed = False
                break
        # Run if not completed (e.g. missing data)
        if completed and len(self._output_paths) > 0:
            logging.info("%s  %s is complete" % (self.record_id, self.name))
        else:
            logging.info("%s  Running %s" % (self.record_id, self.name))
            self.run()

    def _results(self) -> Dict[str, object]:
        # Check that all required datasets are fulfilled
        # Alert if data output is provided, but does not exist
        for _path in self._output_paths.values():
            if isinstance(_path, str) and not os.path.exists(_path):
                # Write dummy file if in developer mode
                if self._mode == 0:
                    touch(_path)
                elif not self.is_skip:
                    print("Missing file: ", _path)
                    raise OutputResultsFileError(_path)
        return self._output_paths

    @property
    def threads(self) -> str:
        return self._threads_pw

    @property
    def program(self) -> LocalCommand:
        if isinstance(self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name], dict):
            return self.local[self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name][ConfigManager.PROGRAM]]
        return self.local[self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name]]

    @property
    def local(self) -> LocalMachine:
        """ Return reference to all commands on user's PATH
        Wrapper for plumbum's LocalMachine object, see plumbum documentation for more info
        https://plumbum.readthedocs.io/en/latest/local_commands.html

        Example: self.local["ls"], self.local["echo"], etc.

        :return: LocalMachine object to use as dictionary to get command to run
        """
        return local

    @property
    def added_flags(self) -> List[str]:
        """ Get additional flags that user provided in config file

        Example: self.local["ls"][(\*self.added_flags("ls"))]

        :return: List of arguments to pass to calling program
        """
        if isinstance(self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name], dict):
            return self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name][ConfigManager.FLAGS].split(" ")
        return []

    @property
    def name(self) -> str:
        """ Name of task and matching config section

        :return: Assigned task name
        """
        return self._name

    @property
    def data(self) -> List[str]:
        return self.config[ConfigManager.DATA].split(",")

    @property
    def output(self) -> Dict[str, Union[Iterable, object]]:
        """ Expected output objects from a successful run
        Output also contains keys that consist of output objects to copy to final results directory

        Example:
        self.output = {"out" : "file-path", "final": ["out"]}
        """
        return self._output_paths

    @output.setter
    def output(self, v: Dict[str, object]):
        self._output_paths = v

    @property
    def input(self) -> Dict[str, Dict[str, Union[Iterable, object]]]:
        """ Dictionary of files available as input to this task.
        By default, available input consists of all files that were passed as input to EukMetaSanity:

        Example:
        self.input["root"] accesses all available input.
        Some available input: self.input["root"]["fna"], self.input["root"]["prot"], self.input["root"]["nr-gff3"]
        This will depend on the pipeline you are using and/or writing

        A task may contain dependencies that must run prior to itself. The output of all of these dependencies
        will be stored within a task's input.

        Example:
        If a task `task2` is dependent on `task`, then the former's input will contain (at least):
        self.input["root"]
        self.input["task"]

        :return: Dictionary of available input data
        """
        return self._input_data

    @property
    def record_id(self) -> str:
        """ Return the ID of a given genome record - typically the basename of the file used in analysis

        :return: ID of file
        """
        return self._record_id

    @property
    def cfg(self) -> ConfigManager:
        """ ConfigManager is a wrapper class for Python's configparser package.
        self.cfg contains references to the entire config file's contents, not just this task's data

        :return: Reference to ConfigManager object that was used to launch this pipeline
        """
        return self._cfg

    @property
    def config(self) -> Dict[str, Union[str, Dict[str, Union[str, dict]]]]:
        """ ConfigManager is a wrapper class for Python's configparser package.
        self.config contains references to the portion of the config file that contains this task's contents

        :return: Dictionary containing contents of config file that was used to launch this task
        """
        if self._scope == "":
            return self._cfg.config[self._name]
        else:
            return self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name]

    @property
    def wdir(self) -> str:
        """ self.wdir contains the directory that should house this task's intermediary files and output
        The directory is automatically created when the task object is created.

        Example:
        self.local["echo"]["Hello world!"] >> os.path.join(self.wdir, "test-output.txt")

        :return: str with the path created for working on current task, storing contents, output, etc.
        """
        return self._wdir

    @property
    def pm(self) -> PathManager:
        """ PathManager object for handling directory creation within task directory

        :return: Reference to PathManager for project
        """
        return self._pm

    def parallel(self, cmd: LocalCommand, time_override: Optional[str] = None):
        """ Launch a command that uses multiple threads
        This method will call a given command on a SLURM cluster automatically (if requested by the user)
        In a config file, WORKERS will correspond to the number of tasks to run in parallel. For slurm users, this
        is the number of jobs that will run simultaneously.

        A time-override may be specified to manually set the maximum time limit a command (job) may run on a cluster,
        which will override the time that is specified by the user in a config file

        The command string will be written to the EukMetaSanity pipeline output file and will be printed to screen

        Example:
        self.parallel(self.local["pwd"], "1:00")

        :param cmd: plumbum LocalCommand object to run
        :param time_override: Time override in "HH:MM:SS" format, if needed
        """
        print("  " + str(cmd))
        # Write command to slurm script file and run
        if self.cfg.config.get(ConfigManager.SLURM, ConfigManager.USE_CLUSTER) is True:
            if ConfigManager.MEMORY not in self.config.keys():
                raise MissingDataError("SLURM section not properly formatted within %s" % self._name)
            cmd = SLURMCaller(
                self.cfg.get_slurm_userid(),
                self.wdir,
                str(self._threads_pw),
                cmd,
                self.config,
                self.local,
                self.cfg.get_slurm_flagged_arguments(),
                time_override
            )
        # Run command directly
        logging.info(str(cmd))
        if self._mode == 1:
            out = cmd()
            if out is not None:
                with open(os.path.join(self.wdir, "task.log"), "a") as w:
                    w.write(str(out))

    def single(self, cmd: LocalCommand):
        """ Launch a command that uses a single thread
        For SLURM users, this method will launch a given command on the node on which EukMetaSanity is launched.

        The command string will be written to the EukMetaSanity pipeline output file and will be printed to screen

        Example:
        self.single(self.local["pwd"])

        :param cmd: plumbum LocalCommand object to run
        """
        print("  " + str(cmd))
        # Run command directly
        logging.info(str(cmd))
        if self._mode == 1:
            out = cmd()
            if out is not None:
                with open(os.path.join(self.wdir, "task.log"), "a") as w:
                    w.write(str(out))

    def create_script(self, cmd: LocalCommand, file_name: str) -> LocalCommand:
        """ Write a command to file and return its value packaged as a LocalCommand.

        This is highly useful when incorporating programs that only launch in the directory in which it was called

        Example:

        script = self.create_script(self.local["ls"]["~"], "cd.sh")

        This will create a file within self.wdir named `cd.sh`, the contents of which will be:

        #!/bin/bash
        cd <wdir> || return

        ls ~


        This can then be run in parallel or singly:
        self.parallel(script)
        self.single(script)

        :param cmd: Command to write to file
        :param file_name: Name of file to create
        """
        _path = os.path.join(self.wdir, file_name)
        fp = open(_path, "w")
        # Write shebang and move to working directory
        fp.write("#!/bin/bash\ncd %s || return\n\n" % self.wdir)
        # Write command to run
        fp.write("".join((str(cmd), "\n")))
        fp.close()
        self.local["chmod"]["+x", _path]()
        return self.local[_path]

    def batch(self, cmds: List[LocalCommand]):
        """ Run a list of commands using self.threads (one task per thread)

        Example:
        self.batch([self.local["pwd"], self.local["echo"]])

        :param cmds: List of LocalCommand objects to run in parallel
        """
        for i in range(0, len(cmds), int(self.threads)):
            running = []
            for j in range(i, i + int(self.threads)):
                if j >= len(cmds):
                    break
                print("  " + str(cmds[j]))
                logging.info(str(cmds[j]))
                if self._mode == 1:
                    f = cmds[j] & BG
                    running.append(f)
            all([_f.wait() for _f in running])


class TaskList(ABC):
    @property
    @abstractmethod
    def requires(cls) -> List[str]:
        pass

    @property
    @abstractmethod
    def depends(cls) -> List[str]:
        pass

    def __init__(self, new_task: type, name: str, cfg: ConfigManager,
                 input_paths: List[Dict[str, Dict[str, object]]],
                 pm: PathManager, record_ids: List[str], mode: int,
                 scope: str):
        # Call data function for pertinent info
        self.name = name
        # Get workers for TaskList
        # Store workers
        self._workers = 1
        self._threads = 1
        self._scope = scope
        if name in cfg.config.keys():
            self._workers = int(cfg.config[name][ConfigManager.WORKERS])
            self._threads = int(cfg.config[name][ConfigManager.THREADS])
        else:
            self._workers = int(cfg.config[scope][ConfigManager.WORKERS])
            self._threads = int(cfg.config[scope][ConfigManager.THREADS])
        # Get log statement
        self._statement = "\nRunning %s protocol using %i worker(s) and %i thread(s) per worker" % (
            self.name, self._workers, self._threads
        )
        # Store list of tasks to complete
        self._tasks: List[Task] = [
            new_task(
                input_path,
                cfg,
                pm,
                record_id,
                name,
                mode,
                scope
            )
            for input_path, record_id in zip(input_paths, record_ids)
        ]
        # Store ConfigManager object
        self._cfg = cfg
        # Store PathManager object
        self._pm = pm
        # Single(0) or Threaded(1) mode
        self._mode = mode

    @property
    def tasks(self) -> List[Task]:
        return self._tasks

    @property
    def scope(self) -> str:
        return self._scope

    def run(self):
        # Single
        logging.info(self._statement)
        for _task in self._tasks:
            if _task.is_skip is False:
                print(self._statement)
                break
        if self._mode == 0:
            for task in self._tasks:
                task._run()
        # Threaded
        else:
            futures = []
            client = Client(n_workers=self._workers, threads_per_worker=1)
            # Run each future
            for _task in self._tasks:
                futures.append(client.submit(_task._run))
            wait(futures)
            client.close()

    def update(self, to_add: List[Dict[str, object]]):
        for task, task_data in zip(self._tasks, to_add):
            task.input.update(task_data)

    def output(self) -> Tuple[List[Dict[str, object]], List[str]]:
        return (
            [task._results() for task in self._tasks],
            [task.record_id for task in self._tasks],
        )


if __name__ == "__main__":
    pass
