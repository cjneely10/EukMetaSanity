"""
Manages the config file, as well as arguments that are set for each part of the pipeline

"""

import os
from pathlib import Path
from typing import List, Tuple, Dict, Union
import yaml
from plumbum import local, CommandNotFound


class InvalidPathError(FileExistsError):
    """ Wraps FileExistsError, raise if user provided a path that does not exist

    """
    pass


class MissingDataError(FileExistsError):
    """ Wraps FileExistsError, raise if user did not fill out a required DATA
    section in a configuration file

    """
    pass


class InvalidProtocolError(ValueError):
    """ Wraps ValueError, raise if user requests to use a protocol that isn't implemented

    """
    pass


class ConfigManager:
    """ ConfigManager handles parsing user-passed config file

    """
    EXPECTED_RESULTS_DIR = "results"
    ROOT = "root"
    SLURM = "SLURM"
    INPUT = "INPUT"
    THREADS = "threads"
    WORKERS = "workers"
    MEMORY = "memory"
    TIME = "time"
    PROTOCOL = "protocol"
    USE_CLUSTER = "USE_CLUSTER"
    DEPENDENCIES = "dependencies"
    PROGRAM = "program"
    FLAGS = "FLAGS"
    DATA = "data"
    BASE = "base"

    def __init__(self, config_path: str):
        """ Create ConfigManager object

        :param config_path: Path to .yaml config file
        """
        self._config = yaml.load(open(str(Path(config_path).resolve()), "r"), Loader=yaml.FullLoader)
        # Confirm all paths in file are valid
        self._validate()

    @property
    def config(self) -> Dict[str, Dict[str, Union[str, dict]]]:
        """ Get config file parsed to dict

        :return: config file parsed to dict
        """
        return self._config

    # Ensure DATA section is valid for all needed databases - mmseqs, etc.
    def _validate(self):
        """ Confirm that data and dependency paths provided in file are all valid.

        Raises MissingDataError

        """
        for task_name, task_dict in self.config.items():
            if "data" in task_dict.keys() and ("skip" not in task_dict.keys() or task_dict["skip"] is not True):
                for _val in task_dict["data"].split(","):
                    if ":" in _val:
                        _val = _val.split(":")[1]
                    if not os.path.exists(str(Path(_val).resolve())):
                        raise MissingDataError("Data for task %s (provided: %s) does not exist!" % (
                            task_name, _val
                        ))
            if "dependencies" in task_dict.keys():
                for prog_name, prog_data in task_dict["dependencies"].items():
                    # Simple - is only a path with no ability to pass flags
                    if isinstance(prog_data, str):
                        try:
                            local[prog_data]
                        except CommandNotFound:
                            # pylint: disable=raise-missing-from
                            raise MissingDataError(
                                "Dependency %s (provided: %s) is not present in your system's path!" % (
                                    prog_name, prog_data))
                    # Provided as dict with program path and FLAGS
                    elif isinstance(prog_data, dict):
                        try:
                            if "program" not in prog_data:
                                # pylint: disable=raise-missing-from
                                raise InvalidPathError(
                                    "Dependency %s is improperly configured in your config file!" % prog_name
                                )
                            bool(prog_data["program"])
                        except CommandNotFound:
                            raise MissingDataError(
                                "Dependency %s (provided: %s) is not present in your system's path!" % (
                                    prog_name, prog_data["program"]))

    def get_slurm_flagged_arguments(self) -> List[Tuple[str, str]]:
        """ Get SLURM arguments from file

        :return: SLURM arguments parsed to input list
        """
        return [
            (key, str(val)) for key, val in self.config["SLURM"].items()
            if key not in {"USE_CLUSTER", "--nodes", "--ntasks", "--mem", "user-id"}
        ]

    def get_slurm_userid(self):
        """ Get user id from slurm section.

        Raises MissingDataError if not present

        :return: user-id provided in config file
        """
        if "user-id" not in self.config["SLURM"].keys():
            raise MissingDataError("SLURM section missing required user data")
        return self.config["SLURM"]["user-id"]
