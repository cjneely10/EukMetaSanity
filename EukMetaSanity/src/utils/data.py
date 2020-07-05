from typing import Dict

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

    # Required for mmseqs taxonomy assignment pipeline
    def taxonomy(self):
        return "taxonomy", self.data["taxonomy"]

    def repeat_modeling(self):
        return "repeats", self.data["repeats"]


if __name__ == "__main__":
    pass
