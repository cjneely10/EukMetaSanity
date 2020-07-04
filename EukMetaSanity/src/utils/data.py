from typing import Dict

"""
Class to handle tracking all required datasets for each task

"""


class Data:
    # Default accessors
    # API
    ACCESS = "-db"
    # Config file
    DATA = "DATA"
    # Input fasta
    IN = "-in"
    # Paths to programs
    PATH = "PATH"
    # Output data
    OUT = "-out"

    def __init__(self):
        self.data: Dict[str, str] = {
            "taxonomy": "ORTHODB",
        }

    # Required for mmseqs taxonomy assignment pipeline
    def taxonomy(self):
        return self.data["taxonomy"]
