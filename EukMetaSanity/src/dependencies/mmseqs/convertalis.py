import os
from typing import List, Union, Type

from yapim import Task, prefix, DependencyInput


class MMSeqsConvertAlis(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "results_files": [
                os.path.join(self.wdir, prefix(database) + "_" + prefix(data) + ".m8")
                for data, database in zip(self.data, self.input["MMSeqsSearch"]["dbs"])]
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("MMSeqsSearch"),
            DependencyInput("MMSeqsCreateDB")
        ]

    def run(self):
        """
        Run mmseqs.convertalis
        """
        # Output results
        for data, database, outfile in zip(self.data, self.input["MMSeqsSearch"]["dbs"],
                                           self.output["results_files"]):
            if not os.path.exists(outfile):
                self.parallel(
                    self.program[
                        "convertalis",
                        str(self.input["MMSeqsCreateDB"]["db"]),  # Input FASTA sequence db
                        data,  # Input augustus-db
                        database,  # Input tax db
                        outfile,  # Output results file
                        "--threads", self.threads,
                        (*self.added_flags)
                    ],
                    "2:00:00"
                )
