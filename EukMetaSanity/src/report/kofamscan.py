from typing import List, Union, Type

from yapim import Task, DependencyInput


class KOFamScan(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "kegg": self.input["KofamscanExecAnnotation"]["kegg"],
            "final": ["kegg"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("KofamscanExecAnnotation")]

    def run(self):
        pass
