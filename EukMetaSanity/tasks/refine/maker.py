import os
from EukMetaSanity import Task, TaskList, program_catch


class MakerIter(TaskList):
    class Maker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        @program_catch
        def run(self) -> None:
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(MakerIter.Maker, "maker", *args, **kwargs)


if __name__ == "__main__":
    pass
