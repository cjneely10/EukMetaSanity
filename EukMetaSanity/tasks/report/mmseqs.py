import os
from EukMetaSanity import Task, TaskList, program_catch


class MMseqsIter(TaskList):
    class MMseqs(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward input
                *[  # BLASTx-like results for each database provided
                    os.path.join(
                        self.wdir, self.record_id + "_results.%s-m8" % os.path.basename(os.path.splitext(db)[0])
                    ) for db in self.data.split(",")
                ]
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            # Generate MMseqs database for proteins
            _file_db = os.path.join(self.wdir, self.record_id + "_db")
            self.log_and_run(
                self.program[
                    "createdb",
                    self.input[0],
                    _file_db,
                ]
            )
            # Search through each database
            for db in self.data.split(","):
                _out_db = _file_db[:-3] + "_%s-results_db" % os.path.basename(os.path.splitext(db)[0])
                # Linear search
                self.log_and_run(
                    self.program[
                        "search",
                        _file_db,  # Input FASTA sequence db
                        db,  # Input augustus-db
                        _out_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        (*self.added_flags),
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
                    ]
                )

    def __init__(self, *args, **kwargs):
        super().__init__(MMseqsIter.MMseqs, "mmseqs", *args, **kwargs)


if __name__ == "__main__":
    pass
