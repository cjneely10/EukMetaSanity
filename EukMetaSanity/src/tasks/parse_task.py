from typing import List
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Helper class will initialize first TaskList object from main function
Outputs match what is required for passing info between Task objects
"""


class Parse:
    def __init__(self, input_files: List[str], cfg: ConfigManager, pm: PathManager,
                 input_prefixes: List[str], debug: int):
        self.input = input_files
        self.cfg = cfg
        self.pm = pm
        self.input_prefixes = input_prefixes
        self.debug = debug
        self.input_files = input_files
        self.input_prefixes = input_prefixes

    def output(self):
        return self
