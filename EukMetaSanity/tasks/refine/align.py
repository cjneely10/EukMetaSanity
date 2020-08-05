import os
from EukMetaSanity import Task, TaskList, program_catch


class AlignIter(TaskList):
    class Align(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".bam"),  # Sorted bam file
            ]

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
                self.program_minimap2[
                    "-ax", self.minimap_option,
                    self.input[4], self.input[-2], self.input[-1],
                ] > os.path.join(self.wdir, self.record_id + ".sam")
            )
            self.log_and_run(
                self.program_sambamba[
                    "view",
                    "-S", os.path.join(self.wdir, self.record_id + ".sam"),
                    "--format=bam", "-o", os.path.join(self.wdir, self.record_id + ".bam"),
                    "-t", self.threads,
                ]
            )
            self.log_and_run(
                self.program_sambamba[
                    "sort",
                    "--format=bam", "-o", os.path.join(self.wdir, self.record_id + ".bam"),
                    "--tmpdir=%s" % os.path.join(self.wdir, "tmp"),
                    "-t", self.threads, "-m", self.memory_limit,
                    os.path.join(self.wdir, self.record_id + ".bam"),
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(AlignIter.Align, "align", *args, **kwargs)


if __name__ == "__main__":
    pass
