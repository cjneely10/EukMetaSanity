import os
from Bio import SeqIO
from EukMetaSanity import Task, TaskList, program_catch


class HHsuiteIter(TaskList):
    class HHsuite(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward prior steps
                *[  # .hhr results for each database provided
                    os.path.join(self.wdir, db, self.record_id + ".hhr")
                    for db in self.data.split(",")
                ]
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            record_fp = SeqIO.parse(self.input[0], "fasta")
            for record in record_fp:
                _file = os.path.join(self.wdir, str(record.id) + ".faa")
                SeqIO.write([record], _file, "fasta")
                # Recruit from uniref
                self.log_and_run(
                    self.program[
                        "-i", _file, "-o", os.path.join(self.wdir, self.record_id + ".uni.hhr"),
                        "-oa3m", os.path.join(self.wdir, self.record_id + ".uni.a3m"),
                        "-cpu", self.threads,
                        "-n", self.numiter,
                        "-d", self.data_uniref,
                        (*self.added_flags),
                    ]
                )
                # Search each database
                for db in self.data.split(","):
                    self.log_and_run(
                        self.program[
                            "-i", os.path.join(self.wdir, self.record_id + ".uni.a3m"),
                            "-o", os.path.join(self.wdir, self.record_id + ".%s.hhr" % db),
                            "-oa3m", os.path.join(self.wdir, self.record_id + ".%s.a3m" % db),
                            "-cpu", self.threads, "-n", "1", "-d", db,
                            (*self.added_flags),
                        ]
                    )

    def __init__(self, *args, **kwargs):
        super().__init__(HHsuiteIter.HHsuite, "hhsuite", *args, **kwargs)


if __name__ == "__main__":
    pass
