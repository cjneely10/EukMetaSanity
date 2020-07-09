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
    @staticmethod
    def _confirm_data(obj):
        # Assert only correct protocols used in config file
        for step, possible_protocols in obj.protocols.items():
            if ConfigManager.PROTOCOL in obj._cfg.config[step].keys():
                assert obj.cfg.config.get(step, ConfigManager.PROTOCOL) in possible_protocols

    def __init__(self, cfg: ConfigManager, name: str):
        self._cfg = cfg
        self._name = name
        # # Update protocols as needed
        self.protocols: Dict[str, Set[str]] = {
            "repeats": {"simple", "full"},
            "abinitio": {"augustus", "gmes"},
        }
        Data._confirm_data(self)

    @property
    def name(self) -> str:
        return self._name

    @property
    def data(self):
        if ConfigManager.DATA in self._cfg.config[self._name].keys():
            return self._cfg.config.get(self._name, ConfigManager.DATA)
        return None

    @property
    def cfg(self) -> ConfigManager:
        return self._cfg

    @added
    def taxonomy(self):
        pass

    @added
    def repeats(self):
        pass

    @added
    def abinitio(self):
        pass

    @added
    def evidence(self):
        pass


if __name__ == "__main__":
    pass
