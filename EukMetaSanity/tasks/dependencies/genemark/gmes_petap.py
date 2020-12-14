import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class GeneMarkPetapIter(TaskList):
    name = "gmes.petap"
    requires = ["taxonomy"]
    depends = ["gmes.prothint"]
    
    class GeneMarkPetap(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                
            }
            
        @program_catch
        def run(self):
            ev_vals = ["--ES"]
            if os.path.getsize(str(self.input["gmes.prothint"]["hints"])) > 0:
                ev_vals = ["--EP", str(self.input["gmes.prothint"]["hints"]),
                           "--evidence", str(self.input["gmes.prothint"]["evidence"])]
            script = self.create_script(
                self.program[
                    "--sequence", str(self.input["root"]["fna"]),
                    (*ev_vals),
                    "--cores", self.threads, (*self.added_flags),
                    ("--fungus" if "fungi" == self.input["taxonomy"]["taxonomy"].kingdom.value.lower() else "")
                ],
                "abinitio.sh"
            )
            # Run script
            self.parallel(script)
            
    def __init__(self, *args, **kwargs):
        super().__init__(GeneMarkPetapIter.GeneMarkPetap, GeneMarkPetapIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
