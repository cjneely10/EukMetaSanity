import os
from EukMetaSanity import Task, TaskList, program_catch


class PfamIter(TaskList):
    class Pfam(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(PfamIter.Pfam, "pfam", *args, **kwargs)


if __name__ == "__main__":
    pass
