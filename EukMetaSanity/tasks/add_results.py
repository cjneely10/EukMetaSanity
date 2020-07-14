import os
from EukMetaSanity import Task, TaskList, program_catch


class AddResultsIter(TaskList):
    class AddResults(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                os.path.join(self.wdir, os.path.basename(_file))
                for _file in self.input
            ]
            print("Evidence")

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            for _file in self.input:
                self.log_and_run(
                    self.local["ln"]["-s", _file, os.path.join(self.wdir, os.path.basename(_file))]
                )

    def __init__(self, *args, **kwargs):
        super().__init__(AddResultsIter.AddResults, "summarize", *args, **kwargs)


if __name__ == "__main__":
    pass
