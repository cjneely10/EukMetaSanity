from typing import List, Union, Type

from yapim import Task, DependencyInput


class AbinitioGeneMark(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "genemark-gff3": self.input["GeneMarkPETAP"]["ab-gff3"],
            "genemark-prot": self.input["GeneMarkPETAP"]["prot"],
            "final": ["genemark-gff3", "genemark-prot"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["Repeats", "Taxonomy"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("GeneMarkPETAP", {"Repeats": {"mask-fna": "fasta"}, "Taxonomy": ["taxonomy"]})
        ]

    def run(self):
        pass
