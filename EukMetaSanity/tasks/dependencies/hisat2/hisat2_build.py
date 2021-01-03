"""
Module holds Hisat2 build functionality
"""

import os
from EukMetaSanity import Task, TaskList
from EukMetaSanity import program_catch, set_complete


class Hisat2BuildIter(TaskList):
    """ TaskList class iterates over Hisat2-build tasks

    name: hisat2.build

    requires: None

    depends: None

    output: db

    expects: fasta

    final: None

    """
    name = "hisat2.build"
    requires = []
    depends = []

    class Hisat2Build(Task):
        """
        Task class calls Hisat2
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "db": os.path.join(self.wdir, self.record_id + "_db")
            }

        @program_catch
        def run(self):
            """
            Run hisat2-build
            """
            self.single(
                self.program[
                    self.dependency_input["fasta"],
                    self.output["db"]
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(Hisat2BuildIter.Hisat2Build, Hisat2BuildIter.name, *args, **kwargs)
