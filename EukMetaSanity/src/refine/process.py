from typing import List, Union, Type

from yapim import Task, DependencyInput


class ProcessMapping(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "sorted_bams": self.input["SambambaSort"]["sorted_bams"]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["MergeSams"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("SambambaSort", {"MergeSams": {"sams": "alignment"}})]

    def run(self):
        pass
