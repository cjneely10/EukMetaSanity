import os
from EukMetaSanity import Task, TaskList, program_catch

"""
Call BRAKER2 pipeline

"""


class BrakerIter(TaskList):
    class Braker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            # Call braker
            self.log_and_run(
                self.program[
                    "--species=%s" % os.path.join(self.wdir, self.record_id),
                    "--genome=%s" % self.input[1],  # Masked genome
                    "--bam=%s" % ",".join()
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(BrakerIter.Braker, "braker", *args, **kwargs)


if __name__ == "__main__":
    pass
