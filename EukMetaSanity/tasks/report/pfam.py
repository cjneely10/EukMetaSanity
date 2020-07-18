import os
from EukMetaSanity import Task, TaskList, program_catch


class PfamIter(TaskList):
    class Pfam(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward initially passed values (in this case prot FASTA file)
                os.path.join(self.wdir, self.record_id + ".tblout"),  # HMM results file
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
                self.program[
                    *self.added_flags,
                    "--cpu", self.threads,
                    "-o", "/dev/null",
                    "--tblout", self.output[-1],
                    self.data, self.input[0],
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(PfamIter.Pfam, "pfam", *args, **kwargs)


if __name__ == "__main__":
    pass
