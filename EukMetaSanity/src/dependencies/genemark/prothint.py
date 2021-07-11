from typing import List, Union, Type

from yapim import Task, DependencyInput, touch


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
        return [
            DependencyInput("MMSeqsFilterTaxSeqDB")
        ]

    def run(self):
        """
        Run gmes.prothint
        """
        try:
            # Run prothint
            self.parallel(
                self.program[
                    str(self.input["fasta"]),
                    str(self.input["MMSeqsFilterTaxSeqDB"]["fastas"][0]),
                    "--workdir", self.wdir,
                    "--threads", self.threads,
                ]
            )
            touch(str(self.output["hints"]))
            touch(str(self.output["evidence"]))
        except:
            touch(str(self.output["hints"]))
            touch(str(self.output["evidence"]))