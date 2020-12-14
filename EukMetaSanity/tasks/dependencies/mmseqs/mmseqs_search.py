import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class SearchIter(TaskList):
    name = "mmseqs.search"
    requires = []
    depends = ["mmseqs.createdb"]
    
    class Search(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(SearchIter.Search, SearchIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
