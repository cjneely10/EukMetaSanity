import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class RepeatMaskerIter(TaskList):
    name = "repmask.repeat_masker"
    requires = ["taxonomy"]
    
    class RepeatMasker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerIter.RepeatMasker, RepeatMaskerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
