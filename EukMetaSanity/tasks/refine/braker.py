import os
from EukMetaSanity import Task, TaskList, program_catch


class BrakerIter(TaskList):
    class Braker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(BrakerIter.Braker, "braker", *args, **kwargs)


if __name__ == "__main__":
    pass
