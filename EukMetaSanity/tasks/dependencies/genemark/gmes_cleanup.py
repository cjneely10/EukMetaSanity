import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class CleanupIter(TaskList):
    name = "gmes.cleanup"
    requires = []
    depends = []
    
    class Cleanup(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(CleanupIter.Cleanup, CleanupIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
