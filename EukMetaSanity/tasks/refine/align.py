import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch

"""
Align RNAseq to genome

"""


class AlignIter(TaskList):
    class Align(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".sorted.bam"),  # Sorted bam file
            ]

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            # Create index
            self.log_and_run(
                self.program_hisat2build[
                    self.input[4],
                    os.path.join(self.wdir, prefix(self.input[4]))
                ]
            )
            # Perform on each transcriptome
            for files_tuple in self.input[-2]:
                _basename = prefix(files_tuple[0])
                # Align to genome
                self.log_and_run(
                    self.program_hisat2[
                        "-p", self.threads,
                        "-x", os.path.join(self.wdir, prefix(self.input[4])),
                        "-1", files_tuple[0],
                        "-2", files_tuple[1],
                        "-S", os.path.join(self.wdir, _basename + ".sam")
                    ]
                )
                # Convert to bam
                self.log_and_run(
                    self.program_sambamba[
                        "view",
                        "-S", os.path.join(self.wdir, _basename + ".sam"),
                        "-f", "bam", "-o", os.path.join(self.wdir, _basename + ".bam"),
                        "-t", self.threads,
                    ]
                )
                # Sort
                self.log_and_run(
                    self.program_sambamba[
                        "sort",
                        "--format=bam", "-o", os.path.join(self.wdir, _basename + ".sorted.bam"),
                        "--tmpdir", os.path.join(self.wdir, "tmp"),
                        "-t", self.threads, "-m", self.memory_limit,
                        os.path.join(self.wdir, _basename + ".bam"),
                    ]
                )
                # Remove intermediaries
                self.local["rm"][os.path.join(self.wdir, _basename + ".bam")]()
                self.local["rm"][os.path.join(self.wdir, _basename + ".sam")]()

    def __init__(self, *args, **kwargs):
        super().__init__(AlignIter.Align, "align", *args, **kwargs)


if __name__ == "__main__":
    pass
