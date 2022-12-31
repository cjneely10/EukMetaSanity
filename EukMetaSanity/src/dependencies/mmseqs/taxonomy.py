import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, Result


class MMSeqsTaxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": os.path.join(self.wdir, f"{self.record_id}.txt"),
            "tax-db": Result(os.path.join(self.wdir, f"{self.record_id}-tax_db"))
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        pass

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsCreateDB")]

    def run(self):
        # Search taxonomy db
        self.parallel(
            self.program[
                "taxonomy",
                self.input["MMSeqsCreateDB"]["db"],
                self.data[0],
                str(self.output["tax-db"]),
                os.path.join(self.wdir, "tmp"),
                (*self.added_flags),
                "--threads", self.threads,
                "--split-memory-limit", str(int(float(self.memory) * 0.7)) + "G",
            ]
        )
        # Generate taxonomy report
        self.single(
            self.program[
                "taxonomyreport",
                self.data[0],
                str(self.output["tax-db"]),
                self.output["tax-report"]
            ],
            "1:00:00"
        )
