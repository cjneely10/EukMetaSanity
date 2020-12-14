import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class GffReadIter(TaskList):
    name = "gmes.gffread"
    requires = []
    depends = ["gmes.gmes_petap"]

    class GffRead(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "gff3": os.path.join(self.wdir, self.record_id + ".gff3")
            }

        @program_catch
        def run(self):
            self.single(
                self.program[
                    self.input["gmes.gmes_petap"]["gtf"],
                    (*self.added_flags),
                    "-o", str(self.output["gff3"])
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(GffReadIter.GffRead, GffReadIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
