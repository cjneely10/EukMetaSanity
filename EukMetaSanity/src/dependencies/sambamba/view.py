import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix


class SambambaView(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "bams": [str(self.wdir.joinpath(prefix(sam_file))) + ".bam" for sam_file in self.input["sams"]]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        for sam_file, bam_file in zip(self.input["sams"], self.output["bams"]):
            if os.path.exists(bam_file):
                continue
            self.parallel(
                self.program[
                    "view",
                    "-S", sam_file,
                    "-t", self.threads,
                    "-o", bam_file,
                    "-f", "bam"
                ]
            )
