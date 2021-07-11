from typing import List, Union, Type

from yapim import Task, DependencyInput


class AbinitioAugustus(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "aug-gff3": self.input["Augustus"]["ab-gff3"],
            "aug-prot": self.input["Augustus"]["prot"],
            "final": ["aug-gff3", "aug-prot"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["Repeats", "Taxonomy"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("Augustus", {"Repeats": ["mask-fna"], "Taxonomy": ["taxonomy"]})
        ]

    def run(self):
        pass
