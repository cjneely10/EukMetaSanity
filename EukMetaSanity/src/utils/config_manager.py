import os
from pathlib import Path
from plumbum import local
from configparser import RawConfigParser
from plumbum.commands.processes import CommandNotFound


class Config(RawConfigParser):
    pass


class InvalidPathError(FileExistsError):
    pass


class MissingDataError(FileExistsError):
    pass


class ConfigManager:
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
        if "ORTHO_DB" not in self.config["DATA"].keys():
            raise MissingDataError("Missing orthodbv10 info!")
        if not os.path.exists(Path(self.config["DATA"]["ORTHO_DB"]).resolve()):
            raise InvalidPathError("Invalid path for orthodbv10")

    # Parse config file for PATH variables and confirm validity
    def _validate_config_paths(self):
        data_in_keys = False
        # Iterate over all values in config file
        for k, value_dict in self.config.items():
            # Ensure all data is valid
            if k == "DATA":
                data_in_keys = True
                self._validate_data()
                continue
            # Ensure PATH sections are valid
            for key in ("PATH",):
                if key in value_dict.keys():
                    try:
                        local[value_dict[key]]()
                    except CommandNotFound:
                        raise InvalidPathError(value_dict[key])
        # Raise error if missing Data section
        if not data_in_keys:
            raise MissingDataError("Missing DATA section!")

    # Gather user-passed flags for analysis
    def get_added_flags(self, _dict_name):
        if "FLAGS" in dict(self.config[_dict_name]).keys():
            return [def_key.lstrip(" ").rstrip(" ")
                    for def_key in self.config[_dict_name]["FLAGS"].rstrip("\r\n").split(",")
                    if def_key != ""]
        else:
            return []


