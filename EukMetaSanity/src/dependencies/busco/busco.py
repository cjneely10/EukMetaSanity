import glob
import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix


class Busco(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "results": self.wdir.joinpath(self.record_id).joinpath("short_summary.txt")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run busco
        """
        results_directory = os.path.dirname(str(self.output["results"]))
        script = self.create_script(
            self.program[
                "-i", self.input["prot"],
                "--out_path", results_directory,
                "-o", prefix(str(self.input["prot"])),
                "-m", self.config["mode"],
                "-l", self.config["lineage"],
                "--cpu", self.threads,
                (*self.added_flags)
            ],
            "busco.sh"
        )
        self.parallel(script)
        # Change name of output file
        os.replace(
            glob.glob(os.path.join(self.wdir, self.record_id, "*", "short_summary*.txt"))[0],
            str(self.output["results"])
        )
