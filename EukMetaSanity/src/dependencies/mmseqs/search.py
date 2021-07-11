import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix


class MMSeqsSearch(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        outfiles = []
        for i, database in enumerate(self.data):
            if database == "":
                continue
            _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, prefix(database) + str(i)))
            outfiles.append(_outfile)
        self.output = {
            "dbs": outfiles
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("MMSeqsCreateDB")
        ]

    def run(self):
        """
        Run mmseqs.search
        """
        for outfile, db_path in zip(self.output["dbs"], self.data):
            if not os.path.exists(outfile + ".index"):
                self.parallel(
                    self.program[
                        self.config["subname"],
                        str(self.input["MMSeqsCreateDB"]["db"]),  # Input FASTA sequence db
                        db_path,  # Input db
                        outfile,  # Output db
                        os.path.join(self.wdir, "tmp"),
                        (*self.added_flags),
                        "--threads", self.threads,
                    ]
                )
