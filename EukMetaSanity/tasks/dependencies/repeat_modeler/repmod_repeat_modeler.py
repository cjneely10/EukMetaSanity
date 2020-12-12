import os
import glob
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound, InvalidPathError, MissingDataError, InvalidProtocolError


class RepeatModelerIter(TaskList):
    name = "repmod.repeat_modeler"
    requires = ["repmod.build_database"]
    
    class RepeatModeler(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            script = self.create_script(
                self.program[
                    "-pa", str(int(self.threads) // 4 or 1),
                    "-engine", "ncbi", "-database", self.input["repmod.build_database"]["db"]
                ],
                "modeler.sh"
            )
            self.parallel(script)
            _output = glob.glob(os.path.join(self.wdir, "RM*"))
            if len(_output) > 0:
                os.replace(
                    _output[0], os.path.join()
                )
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatModelerIter.RepeatModeler, RepeatModelerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
