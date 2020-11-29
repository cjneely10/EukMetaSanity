import os
from EukMetaSanity import Task, TaskList, program_catch

"""
Use search/linsearch in MMseqs program for annotation of all
user-provided datasets

"""


class MMseqsIter(TaskList):
    class MMseqs(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            dbs = (os.path.basename(os.path.splitext(db)[0]) for db in self.data.split(","))
            self.output = [
                *self.input,  # Forward input
                *[  # BLASTx-like results for each database provided
                    os.path.join(
                        self.wdir,
                        self.record_id + "_%s-results.%s-m8" % (db, db)
                    ) for db in dbs
                ]
            ]

        @program_catch
        def run_1(self):
            # Generate MMseqs database for proteins
            _file_db = os.path.join(self.wdir, self.record_id + "_db")
            if not os.path.exists(_file_db):
                self.log_and_run(
                    self.program[
                        "createdb",
                        self.input[0],
                        _file_db,
                    ],
                    "20:00"
                )
            # Search through each database
            for db in self.data.split(","):
                if "rfam" in db:
                    _file_db = _file_db + "_nuc_db"
                    self.log_and_run(
                        self.program[
                            "createdb",
                            self.input[2],
                            _file_db,
                        ],
                        "20:00"
                    )
                search_prog = "search"
                added_flags = self.added_flags
                if os.path.exists(db + ".linidx"):
                    search_prog = "linsearch"
                    for _v in ("-s", "-k"):
                        if _v in added_flags:
                            del added_flags[added_flags.index(_v) + 1]
                            added_flags.remove(_v)
                else:
                    for _v in ("--kmer-per-seq-scale", "--kmer-per-seq"):
                        if _v in added_flags:
                            del added_flags[added_flags.index(_v) + 1]
                            added_flags.remove(_v)
                if "_nuc_db" in _file_db:
                    _file_db = _file_db.replace("_nuc_db", "")
                _out_db = _file_db[:-3] + "_%s-results_db" % os.path.basename(os.path.splitext(db)[0])
                # Linear search
                self.log_and_run(
                    self.program[
                        search_prog,
                        _file_db,  # Input FASTA sequence db
                        db,  # Input augustus-db
                        _out_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        (*added_flags),
                        "--threads", self.threads,
                    ]
                )
                # Output results
                self.log_and_run(
                    self.program[
                        "convertalis",
                        _file_db,  # Input FASTA sequence db
                        db,  # Input augustus-db
                        _out_db,  # Input tax db
                        _out_db[:-3] + ".%s-m8" % os.path.basename(os.path.splitext(db)[0]),  # Output results file
                        "--threads", self.threads,
                    ],
                    "2:00:00"
                )

    def __init__(self, *args, **kwargs):
        super().__init__(MMseqsIter.MMseqs, "mmseqs", *args, **kwargs)


if __name__ == "__main__":
    pass
