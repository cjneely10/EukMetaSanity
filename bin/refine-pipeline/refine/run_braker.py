from typing import List, Union, Type

from yapim import Task, DependencyInput


class RunBraker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "prot": self.input["Braker"]["prot"],
            "gtf": self.input["Braker"]["gtf"],
            "cds": self.input["Braker"]["cds"],
            "final": ["prot", "gtf", "cds"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["ProcessMapping", "GatherProteins"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput(
            "Braker",
            {"ProcessMapping": {"sorted_bams": "bams"}, "GatherProteins": ["prots"], "root": {"mask-fna": "fasta"}}
        )]

    def run(self):
        pass
