import os
from pathlib import Path
from numba import types, jit
from numba.experimental import jitclass
from configparser import RawConfigParser
from configparser import NoSectionError, NoOptionError


class Config(RawConfigParser):
    """ Placeholder class for calling configparser


    """
    pass


class InvalidPathError(FileExistsError):
    pass


class MissingDataError(FileExistsError):
    pass


@jitclass((
    ("_config", Config),
))
class ConfigManager:

    @jit(types.void(types.string), nopython=True, cache=True)
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
    @jit(nopython=True, cache=True)
    def _validate_data(self):
        if "ORTHO_DB" not in self.config.get("DATA", {}).keys():
            raise MissingDataError("Missing orthodbv10 info!")
        if not os.path.exists(Path(self.config["DATA"]["ORTHO_DB"]).resolve()):
            raise InvalidPathError("Invalid path for orthodbv10")

    # Parse config file for PATH variables and confirm validity
    @jit(nopython=True, cache=True)
    def _validate_config_paths(self):
        data_in_keys = False
        # Iterate over all values in config file
        for k, value_dict in self.config.values():
            # Ensure all data is valid
            if k == "DATA":
                data_in_keys = True
                self._validate_data()
                continue
            # Ensure PATH sections are valid
            for key in ("PATH",):
                if key in value_dict.keys() and not os.path.exists(Path(value_dict[key]).resolve()):
                    raise InvalidPathError(value_dict[key])
        # Raise error if missing Data section
        if not data_in_keys:
            raise MissingDataError("Missing DATA section!")

    # Gather user-passed flags for analysis
    @jit(types.ListType(types.string)(types.string), nopython=True, cache=True)
    def get_added_flags(self, _dict):
        if "FLAGS" in dict(self.config[_dict]).keys():
            return [def_key.lstrip(" ").rstrip(" ")
                    for def_key in self.config[_dict]["FLAGS"].rstrip("\r\n").split(",")
                    if def_key != ""]
        else:
            return []


