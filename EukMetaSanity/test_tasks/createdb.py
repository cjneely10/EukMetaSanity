import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class MMSeqsCreateDBIter(TaskList):
    class MMSeqsCreateDB(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(MMSeqsCreateDBIter.MMSeqsCreateDB, "createdb", [],
                         *args, **kwargs)


if __name__ == "__main_":
    pass
