import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class ProtHintIter(TaskList):
    name = "gmes.prothint"
    requires = []
    depends = ["mmseqs.filtertaxseqdb"]
    
    class ProtHint(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(ProtHintIter.ProtHint, ProtHintIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
