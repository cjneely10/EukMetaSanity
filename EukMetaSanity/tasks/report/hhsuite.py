import os
from EukMetaSanity import Task, TaskList, program_catch


class HHsuiteIter(TaskList):
    class HHsuite(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(HHsuiteIter.HHsuite, "hhsuite", *args, **kwargs)


if __name__ == "__main__":
    pass