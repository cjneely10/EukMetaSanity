from typing import List, Union, Type

from yapim import Task, DependencyInput, touch

from .genemark import determine_fungal


class GeneMarkProtHint(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "hints": self.wdir.joinpath("prothint.gff"),
            "evidence": self.wdir.joinpath("evidence.gff")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsFilterTaxSeqDB")]

    def run(self):
        self.parallel(
            self.program[
                str(self.input["fasta"]),
                str(self.input["MMSeqsFilterTaxSeqDB"]["fastas"][0]),
                "--workdir", self.wdir,
                "--threads", self.threads,
                (*determine_fungal(self.input["taxonomy"]))
            ]
        )
        # Create files if prothint fails
        touch(str(self.output["hints"]))
        touch(str(self.output["evidence"]))
