import glob
import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, Result


class GMAPBuild(Task):
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
        Run gmap.build
        """
        _genome_dir = os.path.dirname(str(self.input["fasta"]))
        _genome_basename = os.path.basename(str(self.input["fasta"]))
        if len(glob.glob(str(self.output["db"]) + "*")) == 0:
            self.single(
                self.program[
                    "-d", str(self.output["db"]),
                    "-D", _genome_dir, _genome_basename
                ]
            )
