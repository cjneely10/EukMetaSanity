import os
import glob
import shutil
from EukMetaSanity import Task, TaskList, program_catch, touch


class RepeatModelerIter(TaskList):
    name = "repmod.repeat_modeler"
    requires = ["repmod.build_database"]
    
    class RepeatModeler(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "model": os.path.join(self.wdir, "consensi.fa.classified")
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
            _output = glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa.classified"))
            if len(_output) > 0:
                shutil.move(_output[0], os.path.join(self.wdir, "consensi.fa.classified"))
            else:
                touch(str(self.output["model"]))
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatModelerIter.RepeatModeler, RepeatModelerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
