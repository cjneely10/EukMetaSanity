from typing import Dict, Set

"""
Class to handle tracking all required datasets for each task

"""


class Data:
    # Default accessors
    # API
    ACCESS = "-db"
    # Input fasta
    IN = "-in"
    # Output data
    OUT = "-out"
    # Additional/optional data
    ADDED = "-added"

    def __init__(self):
        self.data: Dict[str, str] = {
            "taxonomy": "ORTHODB",
            "repeats": "MODELER",
        }
        self.protocols: Dict[str, Set[str]] = {
            "repeats": {"simple", "full"},
        }

    # Required for mmseqs taxonomy assignment pipeline
    def taxonomy(self):
        return (
            "taxonomy",
            self.data["taxonomy"],
            "Identifying taxonomy using %i workers and %i threads per worker"
        )

    def repeat_modeling(self):
        return (
            "repeats",
            self.data["repeats"],
            "Identifying repeats using %i workers and %i threads per worker"
        )


if __name__ == "__main__":
    pass
