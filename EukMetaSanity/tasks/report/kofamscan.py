import os
from EukMetaSanity import Task, TaskList, program_catch


class KoFamScanIter(TaskList):
    class KoFamScan(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward initially passed values (in this case prot FASTA file)
                os.path.join(self.wdir, self.record_id + ".kegg.out"),  # KEGG output file
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
                self.program[
                    "-p", self.data_profiles,
                    "-k", self.data_kolist,
                    "--cpu", self.threads,
                    *self.added_flags,
                    "-o", os.path.join(self.wdir, self.record_id + ".kegg.out"),
                    self.input[-1],
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, "kofamscan", *args, **kwargs)


if __name__ == "__main__":
    pass
