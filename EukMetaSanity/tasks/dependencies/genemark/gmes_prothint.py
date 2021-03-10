"""
Module holds gmes.prothint build functionality
"""
import os
from EukMetaSanity import ProcessExecutionError, set_complete
from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput


class ProtHintIter(TaskList):
    """ TaskList class iterates over gmes.prothint tasks

    name: gmes.prothint

    requires:

    depends: mmseqs.filtertaxseqdb

    expects: fasta[Path]

    output: hints[Path], evidence[Path]

    config:
        gmes.prothint:
          program: prothint.py

    """
    name = "gmes.prothint"
    requires = []
    depends = [DependencyInput("mmseqs.filtertaxseqdb")]

    class ProtHint(Task):
        """
        Task class handles gmes.prothint task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "hints": os.path.join(self.wdir, "prothint.gff"),
                "evidence": os.path.join(self.wdir, "evidence.gff")
            }

        @program_catch
        def run(self):
            """
            Run gmes.prothint
            """
            if os.path.getsize(str(self.dependency_input["fasta"])) == 0:
                touch(str(self.output["hints"]))
                touch(str(self.output["evidence"]))
                return
            try:
                # Run prothint
                self.parallel(
                    self.program[
                        str(self.dependency_input["fasta"]),
                        self.input["mmseqs.filtertaxseqdb"]["fastas"][0],
                        "--workdir", self.wdir,
                        "--threads", self.threads,
                    ]
                )
            except ProcessExecutionError:
                touch(str(self.output["hints"]))
                touch(str(self.output["evidence"]))

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(ProtHintIter.ProtHint, ProtHintIter.name, *args, **kwargs)
