from typing import List, Union, Type

from yapim import Task, DependencyInput


class Quality(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "results": self.input["Busco"]["results"],
            "final": ["results"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("Busco")]

    def run(self):
        pass
