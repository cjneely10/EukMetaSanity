import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class FindRNAIter(TaskList):
    class FindRNA(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(FindRNAIter.FindRNA, "rfam", *args, **kwargs)


if __name__ == "__main_":
    pass
