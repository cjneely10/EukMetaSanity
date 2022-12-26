from typing import List, Union, Type

from yapim import Task, DependencyInput


class RunBraker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        output = {}
        final = []
        for key, file in self.input["Braker"]["possible_files"].items():
            if file.exists():
                output[key] = file
                final.append(key)
        self.output = {**output, "final": final}

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
