from typing import List, Union, Type

from yapim import Task, DependencyInput


class RNASeq(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sams": self.input["Hisat2"]["sams"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["CollectInput"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("Hisat2", {"root": {"mask-fna": "fasta"}, "CollectInput": ["rna_read_pairs"]})]

    def run(self):
        pass
