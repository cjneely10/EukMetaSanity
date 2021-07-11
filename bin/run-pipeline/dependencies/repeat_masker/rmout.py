import os
import shutil
from typing import List, Union, Type

from yapim import Task, DependencyInput


class RMaskRMOut(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "mask-gff3": os.path.join(self.wdir, "mask.final.gff3"),
            "mask-fna": os.path.join(self.wdir, self.record_id + ".mask.fna")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("RMaskProcessRepeats")]

    def run(self):
        """
        Run repmask.rmout
        """
        input_file = str(self.input["fasta"]) + ".masked"
        if os.path.exists(input_file):
            os.replace(
                input_file,
                str(self.output["mask-fna"])
            )
        else:
            shutil.copyfile(
                str(self.input["fasta"]),
                str(self.output["mask-fna"])
            )
        # Output the repeats file as a gff3 file
        self.single(
            (self.program[
                 self.input["RMaskProcessRepeats"]["rmout"]
             ] > str(self.output["mask-gff3"]))
        )
