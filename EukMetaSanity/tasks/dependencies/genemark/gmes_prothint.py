"""
Module holds gmes.prothint build functionality
"""
import os
from EukMetaSanity import ProcessExecutionError, set_complete
from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput


class ProtHintIter(TaskList):
    """ TaskList class iterates over gmes.prothint tasks

    name: gmes.prothint

    requires: evidence.prot[Path]

    depends: mmseqs.filtertaxseqdb

    expects: fasta[Path], prot[Path]

    output: hints[Path], evidence[Path]

    config:
        gmes.prothint:
          program: prothint.py

    """
    name = "gmes.prothint"
    requires = ["evidence"]
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
            touch(str(self.output["hints"]))
            touch(str(self.output["evidence"]))
            try:
                tmp_file = os.path.join(self.wdir, "fasta.tmp")
                (self.local["cat"][
                    str(self.input["evidence"]["prot"]), str(self.input["mmseqs.filtertaxseqdb"]["fastas"][0])
                ] > tmp_file)()
                # Run prothint
                self.parallel(
                    self.program[
                        str(self.dependency_input["fasta"]),
                        tmp_file,
                        "--workdir", self.wdir,
                        "--threads", self.threads,
                    ]
                )
            except:
                return

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(ProtHintIter.ProtHint, ProtHintIter.name, *args, **kwargs)
