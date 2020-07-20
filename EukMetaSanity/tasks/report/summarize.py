import os
from EukMetaSanity import Task, TaskList


class SummarizeIter(TaskList):
    class Summarize(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = self.input + ["none"]
            print(self.output)
            for _val in self.output:
                if isinstance(_val, str):
                    print(_val)
                    for _prefix, _ext in os.path.splitext(_val):
                        self.output.append(
                            {_prefix: _ext}
                        )
            print(self.output)

        def run(self):
            super().run()

    def __init__(self, *args, **kwargs):
        super().__init__(SummarizeIter.Summarize, "summarize", *args, **kwargs)


if __name__ == "__main__":
    pass
