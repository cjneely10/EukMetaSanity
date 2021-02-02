"""
Module holds mmseqs.createdb build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class CreateDBIter(TaskList):
    """ TaskList class iterates over mmseqs.createdb tasks

    name: mmseqs.createdb

    requires:

    depends:

    expects: fasta

    output: db

    config:
        mmseqs.createdb:
          program: mmseqs

    """
    name = "mmseqs.createdb"
    requires = []
    depends = []

    class MMSeqsCreateDB(Task):
        """
        Task class handles mmseqs.createdb task
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
            Run mmseqs.createdb
            """
            self.single(
                self.program[
                    "createdb",
                    self.dependency_input["fasta"],
                    self.output["db"]
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(CreateDBIter.MMSeqsCreateDB, CreateDBIter.name, *args, **kwargs)
