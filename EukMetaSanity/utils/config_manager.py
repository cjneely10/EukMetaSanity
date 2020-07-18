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
    # # Protocols for running a choice of a program
    PROTOCOL = "PROTOCOL"
    # # Used for repeating step multiple times
    # ROUNDS = "ROUNDS"

    def __init__(self, config_path):
        self._config = Config()
        self._config.optionxform = str
        self._config.read(config_path)
        # Confirm all paths in file are valid
        self._validate_config_paths()

    @property
    def config(self):
        return self._config

    # Ensure DATA section is valid for all needed databases - mmseqs, etc.
    @staticmethod
    def _validate_data(key, inner_dict):
        for possible_key in (ConfigManager.DATA,):
            _path = inner_dict.get(possible_key, None)
            if _path is not None and not os.path.exists(Path(_path).resolve()):
                raise InvalidPathError("Invalid path %s for %s %s" % (_path, key, possible_key))

    # Parse config file for PATH variables and confirm validity
    def _validate_config_paths(self):
        # Check if protocols are valid
        # Iterate over all values in config file
        for k, value_dict in self.config.items():
            # Ensure all data is valid
            ConfigManager._validate_data(k, value_dict)
            # Ensure PATH sections are valid
            possible_paths = set([val for val in value_dict.keys() if ConfigManager.PROGRAM in val])
            # for possible_path in possible_paths:
            #     if possible_path in value_dict.keys():
            #         try:
            #             local[value_dict[possible_path]]()
            #         except CommandNotFound:
            #             raise InvalidPathError("%s %s" % (k, value_dict[possible_path]))

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


if __name__ == "__main__":
    pass
