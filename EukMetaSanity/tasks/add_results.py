import os
from EukMetaSanity import Task, TaskList, program_catch


class AddResultsIter(TaskList):
    class AddResults(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            for _file in self.input:
                self.log_and_run(
                    self.local["ln"]["-srf", _file, self.wdir]
                )

    def __init__(self, *args, **kwargs):
        super().__init__(AddResultsIter.AddResults, "results", *args, **kwargs)


if __name__ == "__main__":
    pass
