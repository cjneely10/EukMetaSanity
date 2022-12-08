import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix


class GMAP(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sams": [os.path.join(self.wdir, prefix(transcript) + ".sam") for transcript in self.input["transcripts"]]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("GMAPBuild")
        ]

    def run(self):
        """
        Run gmap
        """
        # Get transcripts
        for transcript, sam_file in zip(self.input["transcripts"], self.output["sams"]):
            self.parallel(
                self.program[
                    "-D", str(self.input["GMAPBuild"]["db_dir"]),
                    "-d", str(self.input["GMAPBuild"]["db_name"]),
                    "-t", self.threads,
                    (*self.added_flags),
                    transcript,
                ] > str(sam_file)
            )
