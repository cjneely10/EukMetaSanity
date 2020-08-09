from collections import namedtuple
from typing import Dict, Set, Callable, Tuple
from EukMetaSanity.utils.config_manager import ConfigManager

"""
Class to handle tracking all required datasets for each task

"""

UrlInfo = namedtuple("UrlInfo", ("url", "tar", "flags", "gz"))


def data_urls() -> Dict[str, UrlInfo]:
    return {
        "ortho_db": UrlInfo(
            url="https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz",
            flags="",
            tar=False,
            gz=True,
        )
    }


def added(f: Callable) -> Callable:
    def _wrapper(self) -> Tuple[str, Set[str], str]:
        return (
            f.__name__,
            self.data,
            "\nRunning {} protocol using %i workers and %i threads per worker".format(f.__name__)
        )
    return _wrapper


class Data:
    @staticmethod
    def _confirm_data(obj):
        # Assert only correct protocols used in config file
        for step, possible_protocols in obj.protocols.items():
            if step in obj._cfg.config.keys() and ConfigManager.PROTOCOL in obj._cfg.config[step].keys():
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
        data = []
        for key in self._cfg.config[self._name]:
            if key.startswith(ConfigManager.DATA):
                data.append(self._cfg.config.get(self._name, key))
        return data

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

    @added
    def kofamscan(self):
        pass

    @added
    def mmseqs(self):
        pass

    @added
    def eggnog(self):
        pass

    @added
    def rfam(self):
        pass

    @added
    def summarize(self):
        pass

    @added
    def rnaseq(self):
        pass

    @added
    def transcriptomes(self):
        pass

    @added
    def merge(self):
        pass


if __name__ == "__main__":
    pass
