import os
from EukMetaSanity import Task, TaskList, program_catch


class AssembleIter(TaskList):
    class Assemble(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            files = self.get_reads_files(self.rnaseq)
            built_files = []
            if files is None:
                files = []
            for i in range(0, len(files) - 1, 2):
                built_files.append(
                    os.path.join(
                        self.wdir,
                        self.record_id + "trinity" + os.path.basename(os.path.splitext(files[i])[0]),
                        "Trinity.fasta"
                    )
                )
            self.output = [
                *self.input,
                built_files,  # Assembled transcriptome
                files,  # Identified read files
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            # Get reads files
            files = self.get_reads_files(self.rnaseq)
            for i in range(0, len(files) - 1, 2):
                left, right = files[i], files[i + 1]
                if not (os.path.exists(left) and os.path.exists(right)):
                    continue
                # Run Trinity
                self.log_and_run(
                    self.program[
                        "--seqType", "fq",
                        (*self.added_flags),
                        "--left", left, "--right", right,
                        "--CPU", self.threads,
                        "--output", os.path.join(
                            self.wdir, self.record_id + "trinity" + os.path.basename(os.path.splitext(left)[0])
                        ),
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
