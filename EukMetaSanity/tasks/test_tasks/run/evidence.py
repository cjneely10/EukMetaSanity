import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class EvidenceIter(TaskList):
    name = "evidence"
    requires = []
    
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
