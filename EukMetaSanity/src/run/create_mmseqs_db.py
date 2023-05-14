from typing import List, Union, Type

from yapim import Task, DependencyInput


class CreateMMSeqsDB(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {"db": self.input["MMSeqsCreateDB"]["db"]}

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        pass

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsCreateDB")]

    def run(self):
        pass
