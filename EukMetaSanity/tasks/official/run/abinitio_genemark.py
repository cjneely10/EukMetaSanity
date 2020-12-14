import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class AbInitioGeneMarkIter(TaskList):
    name = "abinitio.genemark"
    requires = ["taxonomy"]
    depends = ["gmes.gffread"]
    
    class AbInitioGeneMark(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "final": "gmes.gffread.gff3"
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioGeneMarkIter.AbInitioGeneMark, AbInitioGeneMarkIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
