from typing import List, Union, Type

from yapim import Task, DependencyInput


class GatherProteins(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "prots": self.input["MMSeqsFilterTaxSeqDB"]["fastas"][0]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["CollectInput"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsFilterTaxSeqDB")]

    def run(self):
        pass
