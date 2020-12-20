import os
from EukMetaSanity import Task, TaskList, program_catch


class GffReadIter(TaskList):
    name = "gmes.gffread"
    requires = []
    depends = ["gmes.petap"]
    
    class GffRead(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": os.path.join(self.wdir, self.record_id + ".gff3")
            }
            
        @program_catch
        def run(self):
            self.single(
                self.program[
                    self.input["gmes.petap"]["gtf"],
                    (*self.added_flags),
                    "-o", str(self.output["ab-gff3"])
                ]
            )
            self.single(
                self.local["sed"][
                    "-i", "s/GeneMark.hmm/ab-initio/g",
                    str(self.output["ab-gff3"])
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(GffReadIter.GffRead, GffReadIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
