import os
from EukMetaSanity import Task, TaskList, program_catch


class MMseqsIter(TaskList):
    """ Task searches gene calls through one (or more) mmseqs databases or profile databases

    Outputs: [db-prefix, dynamic]
    Finalizes: [db-prefix, dynamic]

    """
    name = "mmseqs"
    requires = []
    
    class MMseqs(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            dbs = (os.path.basename(os.path.splitext(db)[0]) for db in self.data.split(","))
            self.output = {
                **{db: os.path.join(
                    self.record_id + "_%s-results.%s-m8" % (db, db)
                ) for db in dbs},
                "final": list(dbs)
            }
            
        @program_catch
        def run(self):
            print(self.input)
            # Generate MMseqs database for proteins
            _file_db = os.path.join(self.wdir, self.record_id + "_db")
            if not os.path.exists(_file_db):
                self.single(
                    self.program[
                        "createdb",
                        self.input["root"]["prot"],
                        _file_db,
                    ]
                )
            # Search through each database
            for db in self.data.split(","):
                if "rfam" in db:
                    _file_db = _file_db + "_nuc_db"
                    self.single(
                        self.program[
                            "createdb",
                            self.input["root"]["fna"],
                            _file_db,
                        ],
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
                self.parallel(
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
                self.parallel(
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
        super().__init__(MMseqsIter.MMseqs, MMseqsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
