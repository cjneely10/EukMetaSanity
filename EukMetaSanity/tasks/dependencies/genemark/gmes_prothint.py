import os
from EukMetaSanity import ProcessExecutionError
from EukMetaSanity import Task, TaskList, program_catch, touch


class ProtHintIter(TaskList):
    name = "gmes.prothint"
    requires = []
    depends = ["mmseqs.filtertaxseqdb"]
    
    class ProtHint(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "hints": os.path.join(self.wdir, "prothint.gff"),
                "evidence": os.path.join(self.wdir, "evidence.gff")
            }
            
        @program_catch
        def run(self):
            try:
                # Run prothint
                self.parallel(
                    self.program[
                        str(self.input["root"]["fna"]),
                        self.input["mmseqs.filtertaxseqdb"]["fasta"],
                        "--workdir", self.wdir,
                        "--threads", self.threads,
                    ]
                )
                os.remove(str(self.input["mmseqs.filtertaxseqdb"]["fasta"]))
            except ProcessExecutionError:
                touch(str(self.output["hints"]))
                touch(str(self.output["evidence"]))
            
    def __init__(self, *args, **kwargs):
        super().__init__(ProtHintIter.ProtHint, ProtHintIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
