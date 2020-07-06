from enum import Enum, auto
from typing import Dict, Set, Callable, Tuple
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Class to handle tracking all required datasets for each task

"""


def added(f: Callable) -> Callable:
    def _wrapper(self) -> Tuple[str, Set[str], str]:
        return (
            f.__name__,
            self.data,
            "Identifying {} using %i workers and %i threads per worker".format(f.__name__)
        )
    return _wrapper


class Data:
    # Default accessors
    # API
    class Type(Enum):
        ACCESS = auto()
        # Input fasta
        IN = auto()
        # Output data
        OUT = auto()
        # Additional/optional data
        ADDED = auto()

    def __init__(self, cfg: ConfigManager, name: str):
        self._cfg = cfg
        self._name = name
        self.protocols: Dict[str, Set[str]] = {
            "repeats": {"simple", "full"},
            "abinitio": {"augustus", "gmes"},
        }
        # Assert only correct protocols used in config file
        for step, possible_protocols in self.protocols.items():
            assert self.cfg.config[step][ConfigManager.PROTOCOL] in possible_protocols

    @property
    def data(self):
        if ConfigManager.DATA in self._cfg.config[self._name]:
            return self._cfg.config.get(self._name, ConfigManager.DATA)
        return None

    @property
    def cfg(self) -> ConfigManager:
        return self._cfg

    # # Add task info below

    @added
    def taxonomy(self):
        pass

    @added
    def repeats(self):
        pass

    @added
    def abinitio(self):
        pass


if __name__ == "__main__":
    pass
