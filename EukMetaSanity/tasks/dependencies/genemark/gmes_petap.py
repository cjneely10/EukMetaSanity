import os
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class GeneMarkPetapIter(TaskList):
    name = "gmes.petap"
    requires = ["taxonomy"]
    depends = [DependencyInput("gmes.prothint")]

    class GeneMarkPetap(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "gtf": os.path.join(self.wdir, "genemark.gtf")
            }

        @program_catch
        def run(self):
            ev_vals = ["--ES"]
            if os.path.getsize(str(self.input["gmes.prothint"]["hints"])) > 0:
                ev_vals = ["--EP", str(self.input["gmes.prothint"]["hints"]),
                           "--evidence", str(self.input["gmes.prothint"]["evidence"])]
            script = self.create_script(
                self.program[
                    "--sequence", str(self.dependency_input["fna"]),
                    (*ev_vals),
                    "--cores", self.threads, (*self.added_flags),
                    ("--fungus"
                     if self.input["taxonomy"]["taxonomy"].kingdom is not None and
                        "fungi" == self.input["taxonomy"]["taxonomy"].kingdom.value.lower() else "")
                ],
                "abinitio.sh"
            )
            # Run script
            self.parallel(script)

    def __init__(self, *args, **kwargs):
        super().__init__(GeneMarkPetapIter.GeneMarkPetap, GeneMarkPetapIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
