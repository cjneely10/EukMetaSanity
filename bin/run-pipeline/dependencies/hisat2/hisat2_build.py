import glob
import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, Result


class Hisat2Build(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "db": Result(os.path.join(self.wdir, self.record_id + "_db"))
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run hisat2-build
        """
        if len(glob.glob(str(self.wdir.joinpath(self.record_id + "_db*")))) > 0:
            return
        self.single(
            self.program[
                self.input["fasta"],
                str(self.output["db"])
            ]
        )
