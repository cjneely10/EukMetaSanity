from typing import List, Union, Type

from yapim import Task, DependencyInput


class MMSeqs(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        final = []
        output = {}
        for i in range(len(self.input["MMSeqsConvertAlis"]["results_files"])):
            final.append(f"results{i}")
            output[f"results{i}"] = self.input["MMSeqsConvertAlis"]["results_files"][i]
        self.output = {
            **output,
            "final": final
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsConvertAlis", {"root": {"prot": "fasta"}})]

    def run(self):
        pass
