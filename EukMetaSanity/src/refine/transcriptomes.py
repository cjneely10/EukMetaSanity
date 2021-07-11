from typing import List, Union, Type

from yapim import Task, DependencyInput


class Transcriptomes(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sams": self.input["GMAP"]["sams"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["CollectInput"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("GMAP", {"root": {"mask-fna": "fasta"}, "CollectInput": ["transcripts"]})]

    def run(self):
        pass
