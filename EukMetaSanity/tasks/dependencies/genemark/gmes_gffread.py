"""
Module holds gmes.gffread build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class GffReadIter(TaskList):
    """ TaskList class iterates over gmes.gffread tasks

    name: gmes.gffread

    requires:

    depends: gmes.petap

    expects:

    output: ab-gff3[Path]

    config:
        gmes.gffread:
          program: gffread
          FLAGS:
            -G

    """
    name = "gmes.gffread"
    requires = []
    depends = [DependencyInput("gmes.petap")]

    class GffRead(Task):
        """
        Task class handles gmes.gffread task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate TaskList
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": os.path.join(self.wdir, self.record_id + ".gmes.gff3")
            }

        @program_catch
        def run(self):
            """
            Run gmes.petap
            """
            self.single(
                self.program[
                    self.input["gmes.petap"]["gtf"],
                    (*self.added_flags),
                    "-o", str(self.output["ab-gff3"])
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(GffReadIter.GffRead, GffReadIter.name, *args, **kwargs)
