from typing import List, Union, Type

from yapim import Task, DependencyInput


class EggNog(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "emapper": self.input["EMapper"]["emapper"],
            "final": ["emapper"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("EMapper")]

    def run(self):
        pass
