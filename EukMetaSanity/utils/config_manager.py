import os
from pathlib import Path
from configparser import RawConfigParser

"""
Manages the config file, as well as arguments that are set for each part of the pipeline

"""


class Config(RawConfigParser):
    pass


class InvalidPathError(FileExistsError):
    pass


class MissingDataError(FileExistsError):
    pass


class ConfigManager:
    # Default accessors
    # Starting PATH value - * * API: THIS IS CALLABLE * *
    PROGRAM = "PROGRAM"
    # Config file
    DATA = "DATA"
    # Workers for task
    WORKERS = "WORKERS"
    # Threads per worker
    THREADS = "THREADS"
    # Memory assigned per worker
    MEMORY = "MEMORY"
    # Time allotted for job
    TIME = "TIME"
    # # Protocols for running a choice of a program
    PROTOCOL = "PROTOCOL"
    # # slurm identifiers
    USE_CLUSTER = "USE_CLUSTER"

    def __init__(self, config_path):
        self._config = Config()
        self._config.optionxform = str
        self._config.read(config_path)
        # Confirm all paths in file are valid
        for k, value_dict in self.config.items():
            # Ensure all data is valid
            ConfigManager._validate_data(k, value_dict)

    @property
    def config(self):
        return self._config

    # Ensure DATA section is valid for all needed databases - mmseqs, etc.
    @staticmethod
    def _validate_data(key, inner_dict):
        for possible_key in (ConfigManager.DATA,):
            _path = inner_dict.get(possible_key, None)
            if _path is not None and not all(
                    [os.path.exists(Path(_p).resolve()) for _p in _path.split(",") if ":" not in _p]):
                raise InvalidPathError("Invalid path %s for %s %s" % (_path, key, possible_key))

    # Gather user-passed flags for analysis
    def get_added_flags(self, _dict_name):
        out = []
        _attrs = set(dir(self))
        for key in dict(self.config[_dict_name]).keys():
            # Parse FLAGS argument from comma-separated
            if key == "FLAGS":
                all(out.append(val) for val in [def_key.lstrip(" ").rstrip(" ")
                                                for def_key in
                                                self.config[_dict_name]["FLAGS"].rstrip("\r\n").split(",")
                                                if def_key != ""])
            # Parse remaining args as dictionary items (for those not used in API)
            # Automatically ignores all capitalized values
            elif key not in dir(self) and not any([key.startswith(_attr) for _attr in dir(self) if _attr.isupper()]) \
                    and not any([key.startswith(_attr) for _attr in self.config[_dict_name] if _attr.isupper()]):
                out.append(key)
                out.append(self.config[_dict_name][key])
        return out

    def get_slurm_flagged_arguments(self):
        return {key: val for key, val in self.config["SLURM"].items()
                # if key not in {ConfigManager.USE_CLUSTER, "FLAGS"}}
                if key not in {ConfigManager.USE_CLUSTER, "--nodes", "--ntasks", "--mem", "user-id"}}

    def get_slurm_userid(self):
        if "user-id" not in self.config["SLURM"]:
            raise MissingDataError("SLURM section missing required user data")
        return self.config["SLURM"]["user-id"]

    def get_FLAGS(self, _dict_name):
        return [val for val in [def_key.lstrip(" ").rstrip(" ")
                                for def_key in
                                self.config[_dict_name]["FLAGS"].rstrip("\r\n").split(",")
                                if def_key != ""]]


if __name__ == "__main__":
    pass
