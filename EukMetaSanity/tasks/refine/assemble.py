import os
from EukMetaSanity import Task, TaskList, program_catch


class AssembleIter(TaskList):
    class Assemble(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            files = self.get_reads_files(self.rnaseq)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + "trinity", "Trinity.fasta"),  # Assembled transcriptome
                *(files if files is not None else [])  # Identified read files
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            # Get reads files
            left, right = self.get_reads_files(self.rnaseq)
            if not (os.path.exists(left) and os.path.exists(right)):
                return
            self.log_and_run(
                self.program[
                    "--seqType", "fq",
                    (*self.added_flags),
                    "--left", left, "--right", right,
                    "--CPU", self.threads,
                    "--output", os.path.join(self.wdir, self.record_id + "trinity"),
                ]
            )

        def get_reads_files(self, file_path: str):
            if not os.path.exists(file_path):
                return "", ""
            with open(file_path, "r") as w:
                for line in w:
                    if line.startswith(self.record_id):
                        return tuple(line.rstrip("\r\n").split("\t")[1].split(","))

    def __init__(self, *args, **kwargs):
        super().__init__(AssembleIter.Assemble, "assemble", *args, **kwargs)


if __name__ == "__main__":
    pass
