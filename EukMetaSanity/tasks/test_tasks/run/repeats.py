import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch


class RepeatsIter(TaskList):
    name = "repeats"
    requires = []
    
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, RepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
