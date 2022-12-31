import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, touch


class RModBuildDatabase(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "db": [os.path.join(self.wdir, self.record_id)],
            "_": self.wdir.joinpath(".done")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run repmod.build_database
        """
        if len(os.listdir(self.wdir)) == 0:
            self.single(
                self.program[
                    "-name", os.path.join(self.wdir, self.record_id),
                    str(self.input["fasta"]),
                ],
                "30:00"
            )
        touch(str(self.output["_"]))
