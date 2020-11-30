import os
from EukMetaSanity import Task, TaskList, program_catch


class EggNogMapperIter(TaskList):
    """ Task runs eggnog-mapper on a set of gene calls

    Outputs: emapper
    Finalizes: emapper

    """
    name = "eggnog"
    requires = []
    
    class EggNogMapper(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "emapper": os.path.join(self.wdir, self.record_id + ".emapper"),
                "final": ["emapper"]
            }
            
        @program_catch
        def run(self):
            self.parallel(
                self.program_python27[
                    self.emapper,
                    "-i", self.input["root"]["prot"],
                    "--output", os.path.join(self.wdir, self.record_id),
                    "--cpu", self.threads,
                    (*self.added_flags),
                ]
            )
            os.replace(
                os.path.join(self.wdir, self.record_id + ".emapper.annotations"),
                str(self.output["emapper"])
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(EggNogMapperIter.EggNogMapper, EggNogMapperIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
