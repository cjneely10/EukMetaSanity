import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class EggNogMapperIter(TaskList):
    name = "eggnog"
    requires = []
    
    class EggNogMapper(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(EggNogMapperIter.EggNogMapper, EggNogMapperIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
