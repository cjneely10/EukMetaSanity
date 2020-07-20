from EukMetaSanity import Task, TaskList


class SummarizeIter(TaskList):
    class Summarize(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = self.input

        def run(self):
            super().run()

    def __init__(self, *args, **kwargs):
        super().__init__(SummarizeIter.Summarize, "summarize", *args, **kwargs)


if __name__ == "__main__":
    pass
