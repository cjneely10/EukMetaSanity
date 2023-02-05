import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, clean


class MMSeqsTaxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": os.path.join(self.wdir, f"{self.record_id}.txt")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        pass

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsCreateDB")]

    @clean("*-tax_db")
    def run(self):
        # Search taxonomy db
        tax_db = os.path.join(self.wdir, f"{self.record_id}-tax_db")
        self.parallel(
            self.program[
                "taxonomy",
                self.input["MMSeqsCreateDB"]["db"],
                self.data[0],
                tax_db,
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
                tax_db,
                self.output["tax-report"]
            ],
            "1:00:00"
        )
