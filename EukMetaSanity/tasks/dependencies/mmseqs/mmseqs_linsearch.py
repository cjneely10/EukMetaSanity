import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class LinSearchIter(TaskList):
    name = "mmseqs.linsearch"
    requires = []
    depends = []
    
    class LinSearch(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(LinSearchIter.LinSearch, LinSearchIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
