import os
from EukMetaSanity import Task, TaskList, program_catch


class KoFamScanIter(TaskList):
    class KoFamScan(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward initially passed values (in this case prot FASTA file)
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
            )

    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, "kofamscan", *args, **kwargs)


if __name__ == "__main__":
    pass
