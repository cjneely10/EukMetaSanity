import glob
import os.path
from typing import List, Union, Type

from yapim import Task, DependencyInput, Result


class GMAPBuild(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        db_name = self.record_id + "_db"
        self.output = {
            "db_name": Result(db_name),
            "db_dir": Result(str(self.wdir)),
            "_": os.path.join(self.wdir, db_name)
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run GMAPBuild
        """
        self.single(
            self.program[
                "-d", str(self.output["db_name"]),
                "-D", str(self.output["db_dir"]),
                str(self.input["fasta"])
            ]
        )
