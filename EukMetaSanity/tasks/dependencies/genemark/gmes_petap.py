import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class GeneMarkPetapIter(TaskList):
    name = "gmes.petap"
    requires = ["taxonomy"]
    depends = []
    
    class GeneMarkPetap(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(GeneMarkPetapIter.GeneMarkPetap, GeneMarkPetapIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
