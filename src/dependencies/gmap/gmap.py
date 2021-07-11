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
            # Generate genome index
            genome_idx = self.input["GMAPBuild"]["db"]
            _genome_dir = os.path.dirname(str(self.input["GMAPBuild"]["db"]))
            _genome_basename = os.path.basename(str(self.input["GMAPBuild"]["db"]))
            # Align
            self.parallel(
                self.program[
                    "-D", _genome_dir, "-d", genome_idx,
                    "-t", self.threads,
                    transcript,
                    (*self.added_flags)
                ] > str(sam_file)
            )
