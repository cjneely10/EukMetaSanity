import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class MergeIter(TaskList):
    class Merge(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(MergeIter.Merge, "merge", *args, **kwargs)


if __name__ == "__main_":
    pass
