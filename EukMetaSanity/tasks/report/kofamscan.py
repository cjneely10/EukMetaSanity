import os
from EukMetaSanity import Task, TaskList, program_catch

"""
Assign KO annotations using kofamscan

"""


class KoFamScanIter(TaskList):
    class KoFamScan(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward initially passed values (in this case prot FASTA file)
                os.path.join(self.wdir, self.record_id + ".kegg"),  # KEGG output file
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
                self.program[
                    "--cpu", self.threads,
                    "--format", "detail",
                    (*self.added_flags),
                    "-o", os.path.join(self.wdir, self.record_id + ".kegg"),
                    "--tmp-dir", os.path.join(self.wdir, "tmp"),
                    self.input[0],
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, "kofamscan", *args, **kwargs)


if __name__ == "__main__":
    pass
