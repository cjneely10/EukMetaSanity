from typing import List, Union, Type

from yapim import Task, DependencyInput


class MergeSams(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        output = self.input["RNASeq"]["sams"]
        output.extend(self.input["Transcriptomes"]["sams"])
        self.output = {
            "sams": output
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["RNASeq", "Transcriptomes"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        pass

    def run(self):
        pass
