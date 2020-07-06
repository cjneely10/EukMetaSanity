import os
from pathlib import Path
from plumbum import local
from configparser import RawConfigParser
from EukMetaSanity.src.utils.data import Data
from plumbum.commands.processes import CommandNotFound

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
    # Paths to programs
    PATH = "PATH"
    PATH2 = "PATH2"
    PATH3 = "PATH3"
    PATH4 = "PATH4"
    # Config file
    DATA = "DATA"
    # Workers for task
    WORKERS = "WORKERS"
    # Threads per worker
    THREADS = "THREADS"
    # Protocols for running a choice of a program
    PROTOCOL = "PROTOCOL"

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
    def _validate_data(self):
        # Taxonomy step
        dt = Data()
        odb = dt.taxonomy()[1]
        if odb not in self.config[ConfigManager.DATA].keys():
            raise MissingDataError("Missing orthodbv10 info!")
        if not os.path.exists(Path(self.config[ConfigManager.DATA][odb]).resolve()):
            raise InvalidPathError("Invalid path for orthodbv10")
        # Repeats modelling step
        modeler_name, modeler_ident = dt.repeats()[0:2]
        if self.config[modeler_name][ConfigManager.PROTOCOL] == "full":
            if modeler_ident not in self.config[ConfigManager.DATA].keys():
                raise MissingDataError("Missing additional repeats data for full search info!")
            path = self.config[ConfigManager.DATA][modeler_ident]
            if path != "None" and not all(os.path.exists(Path(_p).resolve()) for _p in path.split(",") if _p != ""):
                raise InvalidPathError("Invalid paths for additional repeats data")

    # Parse config file for PATH variables and confirm validity
    def _validate_config_paths(self):
        data_in_keys = False
        # Check if protocols are valid
        protocols = Data().protocols
        for step, possible_protocols in protocols.items():
            assert self.config[step][ConfigManager.PROTOCOL] in possible_protocols
        # Iterate over all values in config file
        for k, value_dict in self.config.items():
            # Ensure all data is valid
            if k == ConfigManager.DATA:
                data_in_keys = True
                self._validate_data()
                continue
            # Ensure PATH sections are valid
            for i, possible_path in enumerate((ConfigManager.PATH, ConfigManager.PATH2)):
                if possible_path in value_dict.keys():
                    try:
                        local[value_dict[possible_path]]()
                    except CommandNotFound:
                        raise InvalidPathError(value_dict[possible_path])
        # Raise error if missing Data section
        if not data_in_keys:
            raise MissingDataError("Missing DATA section!")

    # Gather user-passed flags for analysis
    def get_added_flags(self, _dict_name):
        out = []
        for key in dict(self.config[_dict_name]).keys():
            # Parse FLAGS argument from comma-separated
            if key == "FLAGS":
                all(out.append(val) for val in [def_key.lstrip(" ").rstrip(" ")
                                                for def_key in
                                                self.config[_dict_name]["FLAGS"].rstrip("\r\n").split(",")
                                                if def_key != ""])
            # Parse remaining args as dictionary items
            elif key not in (
                ConfigManager.THREADS,
                ConfigManager.WORKERS,
                ConfigManager.DATA,
                ConfigManager.PATH,
                ConfigManager.PATH2,
                ConfigManager.PATH3,
                ConfigManager.PATH4,
            ):
                out.append(key)
                out.append(self.config[_dict_name][key])
        return out


if __name__ == "__main__":
    pass
