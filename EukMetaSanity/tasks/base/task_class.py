"""
Task: Class that manages and handles working directory to complete an operation
TaskList: Collection of Task objects that calls run function on each

"""

import os
import sys
import time
import logging
import traceback
import concurrent.futures
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Callable, Optional, Union, Iterable, Sized, Sequence
# pylint: disable=no-member
from plumbum import colors, local
from plumbum.machines.local import LocalCommand, LocalMachine
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.slurm_caller import SLURMCaller
from EukMetaSanity.tasks.base.dependency_input import DependencyInput
from EukMetaSanity.tasks.base.config_manager import ConfigManager, MissingDataError


def program_catch(func: Callable):
    """ Decorator function checks for ProcessExecutionErrors and FileExistErrors when running an executable
    Logging info recorded and printed to stdout

    :param func: Function to test and whose exceptions to catch
    :return: Decorated function
    """

    def _add_try_except(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        # pylint: disable=broad-except
        except BaseException as err:
            logging.info(err)
            logging.info(traceback.print_exc())
            with open(os.path.join(self.wdir, "task.err"), "a") as w_out:
                w_out.write(str(err) + "\n")
                w_out.write(traceback.format_exc() + "\n")
            print(colors.warn | str(err))

    return _add_try_except


def set_complete(func: Callable):
    """ Decorator function that checks whether a given task is completed. Check on Task object creation, but post-
    child-class self.output update

    :param func: Task class initializer method
    :return: Decorated function, class object modified to store updated self.is_complete status
    """
    def _check_if_complete(self, *args, **kwargs):
        func(self, *args, **kwargs)
        self.is_complete = self.set_is_complete()
    return _check_if_complete


class OutputResultsFileError(FileNotFoundError):
    """ Wrapper class for FileNotFoundError

    """
    pass


# Type alias for allowable input types to task input and dependency_input data members
InputType = Union[object, Iterable, Sized, Sequence]


# pylint: disable=too-many-public-methods
# pylint: disable=invalid-name
class Task(ABC):
    """ Task is an abstract base class that API writers will overwrite to handle Task functionality

    """
    def __init__(self, input_data: Dict[str, Dict[str, InputType]], cfg: ConfigManager, pm: PathManager,
                 record_id: str, db_name: str, mode: int, scope: str,
                 requested_input_data: Dict[str, Dict[str, InputType]],
                 expected_input: Dict[str, InputType]):
        """ Instantiate subclass of Task

        :param input_data: Data dict parsed into self.input
        :param cfg: Reference to global config manager
        :param pm: Reference to global path manager
        :param record_id: Task record id (e.g. basename of input file)
        :param db_name: Tasks scope "_" name
        :param mode: Run/debug mode, set at command line level
        :param scope: Outer scope. If not present, Task is outer-level task. Else, is a dependency
        :param requested_input_data: Additional data requested from Task depends and requires lists
        :param expected_input: Input data passed as dependency_input to Task dependencies
        """
        # Set name of task
        self._name = db_name
        # Store input data and update passed/requested data
        self._input_data = input_data
        self._input_data.update(requested_input_data)
        # Store override input for dependency
        self._dep_input = expected_input
        # Instantiate output dict variable
        self._output_paths: Dict[str, object] = {}
        # Store config manager
        self._cfg = cfg
        # Store threads and scope
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
        # Create working directory based on name of task and scope
        db_name = scope + "_" + db_name if scope != "" else db_name
        pm.add_dirs(record_id, [db_name])
        # Store working directory
        self._wdir = pm.get_dir(record_id, db_name)
        # Store id of record in Task
        self._record_id = record_id
        # Check if task is set to be skipped in config file
        # Outer level task - set based on value
        if self._scope == "":
            self.is_skip = "skip" in self.config.keys() and self.config["skip"] is True
        else:
            # Inner level task may have skip set in its config section
            if "skip" in self.config.keys():
                self.is_skip = self.config["skip"] is True
            # Otherwise base on skip status of parent dependency
            else:
                self.is_skip = "skip" in self._cfg.config[self._scope].keys() and \
                           self._cfg.config[self._scope]["skip"] is True
        # Store whether this task has been completed
        self.is_complete = False
        super().__init__()

    @abstractmethod
    def run(self):
        """ Abstract method to run for each task.

        """
        pass

    def run_helper(self) -> None:
        """ Helper method to determine if Task needs is scheduled to run (e.g. not is_skip)
        and is not already completed. Prints completion status to stdout

        """
        # If task is set to skip,return
        if self.is_skip:
            return
        if self.is_complete:
            _str = "Is complete: {}".format(self.record_id)
            logging.info(_str)
            print(colors.blue & colors.bold | _str)
        else:
            _str = "In progress:  {}".format(self.record_id)
            logging.info(_str)
            print(colors.blue & colors.bold | _str)
            start_time = time.time()
            self.run()
            end_time = time.time()
            _str = "Is complete:  {} ({:.3f}{})".format(self.record_id, *Task._parse_time(end_time - start_time))
            logging.info(_str)
            print(colors.blue & colors.bold | _str)

    def is_non_file_output(self) -> bool:
        """ Checks if task outputs and output files at all or is simply some wrapper task

        :return: True if no string value is present in task's output
        """
        for _path in self._output_paths.values():
            if isinstance(_path, str):
                return False
        return True

    def set_is_complete(self) -> bool:
        """ Check all required output data to see if any part of task need to be completed

        :return: Boolean representing if task has all required output
        """
        is_complete = None
        for _path in self._output_paths.values():
            if isinstance(_path, str):
                if not os.path.exists(_path):
                    # Only call function if missing path
                    # Then move on
                    is_complete = False
                    break
                is_complete = True
        if is_complete is None:
            return False
        return is_complete

    def results(self) -> Dict[str, InputType]:
        """ Check that all required output is created

        :raises: OutputResultsFileError if a file is expected to have generated but didn't
        :return: Results dictionary that was set in Task __init__
        """
        for _path in self._output_paths.values():
            if isinstance(_path, str) and not os.path.exists(_path):
                # Write dummy file if in developer mode
                if self._mode == 0:
                    continue
                if not self.is_skip:
                    print(colors.red & colors.bold | "Missing file: %s" % _path)
                    raise OutputResultsFileError(_path)
        return self._output_paths

    @property
    def threads(self) -> str:
        """ Number of threads when running task (as set in config file)

        :return: Str of number of tasks
        """
        return self._threads_pw

    @property
    def program(self) -> LocalCommand:
        """ Run program assigned to task. This will be attached to a given config file section
        based on its outer scope

        :return: LocalCommand to `program` path
        """
        return self.local[self.config[ConfigManager.PROGRAM]]

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

        Example: self.local["ls"][(*self.added_flags("ls"))]

        :return: List of arguments to pass to calling program
        """
        if isinstance(self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name], dict):
            if ConfigManager.FLAGS in self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name].keys():
                out = self._cfg.config[self._scope][ConfigManager.DEPENDENCIES][self._name][
                    ConfigManager.FLAGS].split(" ")
                while "" in out:
                    out.remove("")
                return out
        return []

    @property
    def name(self) -> str:
        """ Name of task and matching config section

        :return: Assigned task name
        """
        return self._name

    @property
    def data(self) -> List[str]:
        """ List of data files passed to this task's config section

        :return: List of paths of data in task's config section
        """
        return self.config[ConfigManager.DATA].split(" ")

    @property
    def output(self) -> Dict[str, InputType]:
        """ Expected output objects from a successful run
        Output also contains keys that consist of output objects to copy to final results directory

        Example:
        self.output = {"out" : "file-path", "final": ["out"]}
        """
        return self._output_paths

    @output.setter
    def output(self, v: Dict[str, InputType]):
        """ Dict of data that is output by this task

        :param v: Dict of str: object that will output when the task successfully completes
        """
        self._output_paths = v

    @property
    def dependency_input(self) -> Dict[str, InputType]:
        """ Input to a dependency. Used to run a dependency using the output of a separate abstract
        Task output

        :return: Dict of str:object of the data that a dependency is expecting as input
        """
        return self._dep_input

    @property
    def input(self) -> Dict[str, Dict[str, InputType]]:
        """ Dictionary of files available as input to this task.
        By default, available input consists of all files that were passed as input to EukMetaSanity:

        Example:
        self.input["root"] accesses all available input.
        Some available input: self.input["root"]["fasta"], self.input["root"]["prot"], self.input["root"]["nr-gff3"]
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
        """ Task's corresponding config file section as a dictionary

        :return: Dictionary containing contents of config file that was used to launch this task
        """
        # Task is an outer-level abstract task
        if self._scope == "":
            return self._cfg.config[self._name]
        # Task is an inner-level dependency of an abstract task
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

    @property
    def developer_mode(self) -> bool:
        """ Return if running in developer mode - only run through metadata, do not run actual methods

        :return: True if in developer_mode
        """
        return bool(self._mode) is False

    @developer_mode.setter
    def developer_mode(self, mode: bool):
        """ Set developer mode (True) or regular run mode (False

        :param mode: True to set developer mode, False otherwise
        """
        # 0 is developer, 1 is user
        # True is False = 0 for developer
        # False is False = 1 for user
        self._mode = int(mode is False)

    @property
    def memory(self) -> str:
        """ Get SLURM memory set at either the current task level, or its parent scope

        :return: Memory string
        """
        if ConfigManager.MEMORY in self.config.keys():
            return self.config[ConfigManager.MEMORY]
        return self.cfg.config[self._scope][ConfigManager.MEMORY]

    @property
    def time(self) -> str:
        """ Get SLURM time set at either the current task level, or its parent scope

        :return: Time string
        """
        if ConfigManager.TIME in self.config.keys():
            return self.config[ConfigManager.TIME]
        return self.cfg.config[self._scope][ConfigManager.TIME]

    @property
    def is_slurm(self) -> bool:
        """
        Run was launched on SLURM
        """
        return bool(self.cfg.config.get(ConfigManager.SLURM)[ConfigManager.USE_CLUSTER])

    def _create_slurm_command(self, cmd: LocalCommand, time_override: Optional[str] = None,
                              threads_override: str = None, memory_override: str = None) -> SLURMCaller:
        """ Create a SLURM-managed process

        :param cmd: plumbum LocalCommand object to run
        :param time_override: Time override in "HH:MM:SS" format, if needed
        :param threads_override: Provide number of threads to parallelize over, default to use config-level threads-pw.
            Note that this will only affect SLURM script generation - this will not override thread values passed in by
            the cmd parameter
        :param memory_override: Provide memory override for command in "2GB" format, etc.
        :return: SLURM-wrapped command to run script via plumbum interface
        """
        sel = self._scope if (self._scope is not None and self._scope != "") else self._name
        # Confirm valid SLURM section
        if ConfigManager.MEMORY not in self.cfg.config[sel].keys():
            raise MissingDataError("SLURM section not properly formatted within %s" % self._name)
        if ConfigManager.TIME not in self.cfg.config[sel].keys():
            raise MissingDataError("SLURM section not properly formatted within %s" % self._name)
        # Generate command to launch SLURM job
        return SLURMCaller(
            self.cfg.get_slurm_userid(),
            self.wdir,
            str(self._threads_pw) if threads_override is None else threads_override,
            cmd,
            self.memory if memory_override is None else memory_override,
            self.time if time_override is None else time_override,
            self.local,
            self.cfg.get_slurm_flagged_arguments(),
        )

    def parallel(self, cmd: Union[LocalCommand, List[LocalCommand]], time_override: Optional[str] = None,
                 threads_override: str = None, memory_override: str = None):
        """ Launch a command that uses multiple threads
        This method will call a given command on a SLURM cluster automatically (if requested by the user)
        In a config file, WORKERS will correspond to the number of tasks to run in parallel. For slurm users, this
        is the number of jobs that will run simultaneously.

        A time-override may be specified to manually set the maximum time limit a command (job) may run on a cluster,
        which will override the time that is specified by the user in a config file

        The command string will be written to the EukMetaSanity pipeline output file and will be printed to screen

        Example:
        self.parallel(self.local["pwd"], "1:00")

        :param cmd: plumbum LocalCommand object to run, or list of commands to run
        :param time_override: Time override in "HH:MM:SS" format, if needed
        :param threads_override: Provide number of threads to parallelize over, default to use config-level threads-pw.
            Note that this will only affect SLURM script generation - this will not override thread values passed in by
            the cmd parameter
        :param memory_override: Provide memory override for command in "2GB" format, etc.
        :raises: MissingDataError if SLURM section improperly configured
        """
        # Write command to slurm script file and run
        if self.is_slurm:
            cmd = self._create_slurm_command(cmd, time_override, threads_override, memory_override)
        # Run command directly
        if self._mode == 1:
            logging.info(str(cmd))
            print("  " + str(cmd))
            out = cmd()
            # Store log info in any was generated
            if out is not None:
                with open(os.path.join(self.wdir, "task.log"), "a") as w:
                    w.write(str(out))

    def single(self, cmd: Union[LocalCommand, List[LocalCommand]],
               time_override: Optional[str] = None, memory_override: str = None):
        """ Launch a command that uses a single thread.

        The command string will be written to the EukMetaSanity pipeline output file and will be printed to screen

        Example:
        self.single(self.local["pwd"])

        :param cmd: plumbum LocalCommand object to run, or list of commands to run
        :param time_override: Time override in "HH:MM:SS" format, if needed
        :param memory_override: Provide memory override for command in "2GB" format, etc.
        """
        self.parallel(cmd, time_override, threads_override="1", memory_override=memory_override)

    def create_script(self, cmd: Union[str, LocalCommand, List[str]], file_name: str) -> LocalCommand:
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

        :param cmd: Command to write to file, or list of commands to write
        :param file_name: Name of file to create
        :return: Command to run script via plumbum interface
        """
        _path = os.path.join(self.wdir, file_name)
        fp = open(_path, "w")
        # Write shebang and move to working directory
        fp.write("#!/bin/bash\ncd %s || return\n\n" % self.wdir)
        # Write command to run
        if isinstance(cmd, list):
            for _cmd in cmd:
                fp.write("".join((str(_cmd), "\n")))
        else:
            fp.write("".join((str(cmd), "\n")))
        fp.close()
        self.local["chmod"]["+x", _path]()
        return self.local[_path]

    def batch(self, cmds: List[LocalCommand], time_override: Optional[str] = None):
        """ Run a list of commands using self.threads (one task per thread)

        Example:
        self.batch([self.local["pwd"], self.local["echo"]])

        :param cmds: List of LocalCommand objects to run in parallel
        :param time_override: Time override in "HH:MM:SS" format, if needed
        """
        write_string = f"  Running {len(cmds)} commands similar to {str(cmds[0])}"
        print(write_string)
        logging.info(write_string)
        if self.is_slurm:
            for i in range(0, len(cmds), int(self.threads)):
                running = []
                # Run up to `threads` tasks at a time
                for j in range(i, i + int(self.threads)):
                    if j >= len(cmds):
                        break
                    running.append(cmds[j])
                self.parallel(running, time_override=time_override)
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=int(self.threads)) as executor:
                running = []
                for cmd in cmds:
                    running.append(executor.submit(self.single, cmd, time_override))
                concurrent.futures.wait(running)

    @staticmethod
    def _parse_time(_time: float) -> Tuple[float, str]:
        """ Parse time to complete a task into
        day, hour, minute, or second representation based on scale

        :param _time: Time to complete a given task
        :return: time and string representing unit
        """
        if _time > 3600 * 24:
            return _time / (3600 * 24), "d"
        if _time > 3600:
            return _time / 3600, "h"
        if _time > 60:
            return _time / 60, "m"
        return _time, "s"


class TaskList(ABC):
    """ TaskList is an abstract base class that API writers will define. The class handles distributing Task and
    ensuring that Tasks are only run as needed, and that Dask is not called to unnecessarily distribute
    Tasks

    """
    @property
    @abstractmethod
    def requires(self) -> List[str]:
        """ Requirements to run a given TaskList

        :return: List of task class names that must run before this task
        """
        pass

    @property
    @abstractmethod
    def depends(self) -> List[DependencyInput]:
        """ Dependencies used to run task

        :return: List of dependencies and the task output to use as its input.
        This defaults to `root` info passed from the pipeline's initialization unless
        explicitly specified in DependencyInput constructor
        """
        pass

    def __init__(self, new_task: type, name: str, cfg: ConfigManager,
                 input_paths: List[Dict[str, Dict[str, InputType]]],
                 pm: PathManager, record_ids: List[str], mode: int,
                 scope: str, requested_input_data: List[Dict[str, Dict[str, InputType]]],
                 expected_input_list: List[Dict[str, InputType]]):
        """ Instantiate child class of TaskList with provided Task type and name. Pass additional input
        requested at Task Level

        :param new_task: Task type for each task in this TaskList collection
        :param name: Name of Task as specified in config file
        :param cfg: Reference to ConfigManager
        :param input_paths: List of input data/path data needed for each Task
        :param pm: Reference to PathManager
        :param record_ids: List of record_ids associated with stored list of tasks
        :param mode: Run/debug mode
        :param scope: Outer scope of Task, may be "" if Task is at outer scope
        :param requested_input_data: List of requested input data for Task that is fed to each Task
        :param expected_input_list: List of dependency input data needed for dependency-level input
        """
        # Name of TaskList
        self.name = name
        # Store workers, threads, and scope
        self._workers = 1
        self._threads = 1
        self._scope = scope
        self.full_name = self.scope + ("_%s" % self.name if self.scope != "" else self.name)
        if name in cfg.config.keys():
            self._workers = int(cfg.config[name][ConfigManager.WORKERS])
            self._threads = int(cfg.config[name][ConfigManager.THREADS])
        else:
            self._workers = int(cfg.config[scope][ConfigManager.WORKERS])
            self._threads = int(cfg.config[scope][ConfigManager.THREADS])
        # Get log statement
        self._statement = "\nRunning:\n  %s\n  %i worker(s)\n  %i thread(s) per worker" % (
            self.full_name, self._workers, self._threads
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
                scope,
                req_data,
                exp_input
            )
            for input_path, record_id, req_data, exp_input in zip(
                input_paths, record_ids, requested_input_data, expected_input_list
            )
        ]
        # Store ConfigManager object
        self._cfg = cfg
        # Store PathManager object
        self._pm = pm
        # Single(0) or Threaded(1) mode
        self._mode = mode

    @property
    def tasks(self) -> List[Task]:
        """ List of Task objects that this TaskList manages

        :return: List of Tasks managed by this TaskList object
        """
        return self._tasks

    @property
    def scope(self) -> str:
        """ Scope assigned to tasklist. For abstract tasks, the scope is "".
        For dependencies, the scope is the name of the abstract task that called it.

        :return: Name of scope, "" if an outer Task, or an outer Task name if is an inner dependency
        """
        return self._scope

    def run(self):
        """ Run TaskList on all Tasks that have not been completed. Distribute over SLURM/parallel
        as needed

        """
        # Confirm that task is not meant to be skipped
        is_skip = True
        for _task in self._tasks:
            if _task.is_skip is False:
                print(colors.green & colors.bold | self._statement)
                is_skip = False
                break
        if is_skip is True:
            return
        # Check if all tasks have been completed
        if self._prerun_check() is True:
            return
        # Log task beginning
        logging.info(self._statement)
        # Run each task in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=self._workers) as executor:
            output_futures = [executor.submit(_task.run_helper) for _task in self._tasks]
            concurrent.futures.wait(output_futures)

    def _prerun_check(self) -> bool:
        """ Check if each task has been completed.

        :return: True if all Tasks in TaskList done, else False
        """
        for _task in self._tasks:
            if not _task.is_complete:
                return False
        return True

    def output(self) -> Tuple[List[Dict[str, object]], List[str]]:
        """ Get output and record_id that is generated by all Tasks in TaskList

        :return: List of output from all Tasks in TaskList and their accompanying record_ids
        """
        try:
            return [task.results() for task in self._tasks], [task.record_id for task in self._tasks]
        except OutputResultsFileError:
            sys.exit()
