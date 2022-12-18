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
                # Profile databases should not use linsearch
                if "p:" in db_path:
                    db_path = db_path[2:]
                    subname = "search"
                else:
                    subname = self.config["subname"]
                self.parallel(
                    self.program[
                        subname,
                        str(self.input["MMSeqsCreateDB"]["db"]),  # Input FASTA sequence db
                        db_path,  # Search db
                        outfile,  # Output db
                        os.path.join(self.wdir, "tmp"),
                        "--split-memory-limit", str(int(float(self.memory) * 0.7)) + "G",
                        (*self.added_flags),
                        "--threads", self.threads,
                    ]
                )
