from enum import Enum, auto
from typing import Dict, Set, Callable

"""
Class to handle tracking all required datasets for each task

"""


def added(f: Callable):
    def _wrapper(self):
        return (
            f.__name__,
            (self.data[f.__name__] if f.__name__ in self.data.keys() else None),
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

    def __init__(self):
        self.data: Dict[str, str] = {
            "taxonomy": "ORTHODB",
            "repeats": "MODELER",
        }
        self.protocols: Dict[str, Set[str]] = {
            "repeats": {"simple", "full"},
            "abinitio": {"augustus", "gmes"},
        }

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
