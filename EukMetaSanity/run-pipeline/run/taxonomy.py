from typing import List, Union, Type

from yapim import Task, DependencyInput


class Taxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "taxonomy": self.input["MMSeqsTaxonomy"]["taxonomy"],
            "taxonomy-actual": self.input["MMSeqsTaxonomy"]["taxonomy-actual"],
            "tax-report": self.input["MMSeqsTaxonomy"]["tax-report"],
            "final": ["tax-report", "taxonomy"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["Evidence"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("MMSeqsTaxonomy", {"Evidence": {"evidence-prot": "fasta"}})
        ]

    def run(self):
        """
        Run
        """
        pass
