import os
from EukMetaSanity import Task, TaskList, program_catch


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
            # Align to genome
            for i in range(0, len(self.input[-1]), 2):
                self.log_and_run(
                    self.program_minimap2[
                        "-ax", self.minimap_option,
                        self.input[4], self.input[-1][i], self.input[-1][i + 1],
                        (*self.added_flags),
                    ] > os.path.join(self.wdir, self.record_id + "%s.sam" % str(i))
                )
                # Convert to bam
                self.log_and_run(
                    self.program_sambamba[
                        "view",
                        "-S", os.path.join(self.wdir, self.record_id + "%s.sam" % str(i)),
                        "--format=bam", "-o", os.path.join(self.wdir, self.record_id + "%s.bam" % str(i)),
                        "-t", self.threads,
                    ]
                )
                # Sort
                self.log_and_run(
                    self.program_sambamba[
                        "sort",
                        "--format=bam", "-o", os.path.join(self.wdir, self.record_id + "%s.sorted.bam" % str(i)),
                        "--tmpdir=%s" % os.path.join(self.wdir, "tmp"),
                        "-t", self.threads, "-m", self.memory_limit,
                        os.path.join(self.wdir, self.record_id + "%s.bam" % str(i)),
                    ]
                )
                # Remove intermediaries
                self.local["rm"][os.path.join(self.wdir, self.record_id + "%s.bam" % str(i))]()
                self.local["rm"][os.path.join(self.wdir, self.record_id + "%s.sam" % str(i))]()

    def __init__(self, *args, **kwargs):
        super().__init__(AlignIter.Align, "align", *args, **kwargs)


if __name__ == "__main__":
    pass
