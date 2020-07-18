import os
from EukMetaSanity import TaskList, Task, program_catch


class EggNOGMapper(TaskList):
    class EggNOG(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,  # Forward values
                os.path.join(self.wdir, self.record_id + ".emapper.annotations")  # Results of mapper
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.log_and_run(
                self.program_python27[
                    self.emapper,
                    "-i", self.input[0],
                    "--output", self.output[-1],
                    "--cpu", self.threads,
                    (*self.added_flags),
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(EggNOGMapper.EggNOG, "eggnog", *args, **kwargs)


if __name__ == "__main__":
    pass
