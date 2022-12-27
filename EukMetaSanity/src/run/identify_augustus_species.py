import os
from typing import List, Union, Type

from yapim import Task, DependencyInput


class IdentifyAugustusSpecies(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "search_results": self.input["MMSeqsConvertAlis"]["results_files"][0],
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["AbinitioGeneMark"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsConvertAlis")]

    def condition(self) -> bool:
        genemark_output = str(self.input["AbinitioGeneMark"]["genemark-gff3"])
        return (not os.path.exists(genemark_output)) or os.stat(genemark_output).st_size == 0

    def run(self):
        pass
