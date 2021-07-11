import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix

try:
    from sambamba.utils import is_sam
except ImportError:
    from tests.eukmetasanity.dependencies.sambamba.utils import is_sam


class SambambaSort(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sorted_bams": [self.wdir.joinpath(prefix(bam_file) + ".sorted.bam")
                            for bam_file in self.input["alignment"]]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run sambamba.sort
        """
        for alignment_file, sorted_bam_file in zip(self.input["alignment"], self.output["sorted_bams"]):
            if os.path.exists(sorted_bam_file):
                continue
            if is_sam(alignment_file):
                self.from_sam(str(alignment_file), str(sorted_bam_file))
            else:
                self.from_bam(str(alignment_file), str(sorted_bam_file))

    def from_sam(self, sam_file: str, sorted_bam_file: str):
        self.parallel(
            (self.program[
                "view",
                "-S", sam_file,
                "-t", self.threads,
                "-o", "/dev/stdout",
                "-f", "bam"
            ] |
             self.program[
                 "sort",
                 "-t", self.threads,
                 "-o", sorted_bam_file,
                 "-m", str(self.memory) + "GB",
                 "/dev/stdin",
                 (*self.added_flags)
             ])
        )

    def from_bam(self, bam_file: str, sorted_bam_file: str):
        self.parallel(
            self.program[
                "sort",
                "-t", self.threads,
                "-o", sorted_bam_file,
                "-m", str(self.memory) + "GB",
                bam_file,
                (*self.added_flags)
            ]
        )
