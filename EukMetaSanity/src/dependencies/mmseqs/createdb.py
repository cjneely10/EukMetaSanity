import os
from typing import List, Union, Type

from yapim import Task, DependencyInput


class MMSeqsCreateDB(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "db": os.path.join(self.wdir, self.record_id + "_db")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run mmseqs.createdb
        """
        self.single(
            self.program[
                "createdb",
                self.input["fasta"],
                self.output["db"]
            ],
            "30:00",
        )
