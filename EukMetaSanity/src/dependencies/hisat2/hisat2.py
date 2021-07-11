import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix


class Hisat2(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sams": [os.path.join(self.wdir, prefix(pair[0]) + ".sam") for pair in self.input["rna_read_pairs"]]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("Hisat2Build")
        ]

    def run(self):
        """
        Run hisat2
        """
        for pair, out_sam in zip(self.input["rna_read_pairs"], self.output["sams"]):
            if os.path.exists(out_sam):
                continue
            if len(pair) > 1:
                _parse_args = ["-1", pair[0], "-2", pair[1]]
            else:
                _parse_args = ["-U", pair[0]]
            # Align
            self.parallel(
                self.program[
                    "-p", self.threads,
                    "-x", self.input["Hisat2Build"]["db"],
                    (*_parse_args),
                    "-S", out_sam,
                    (*self.added_flags),
                ]
            )
