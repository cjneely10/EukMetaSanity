"""
Module holds repmod.build_database build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class BuildDatabaseIter(TaskList):
    """ TaskList class iterates over repmod.build_database tasks

    name: repmod.build_database

    requires:

    depends:

    expects: fasta[Path]

    output: db[Path]

    config:
        repmod.build_database:
          program: BuildDatabase

    """
    name = "repmod.build_database"
    requires = []
    depends = []

    class BuildDatabase(Task):
        """
        Task class handles repmod.build_database task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "db": [os.path.join(self.wdir, self.record_id)]
            }

        @program_catch
        def run(self):
            """
            Run repmod.build_database
            """
            if len(os.listdir(self.wdir)) == 0:
                self.single(
                    self.program[
                        "-name", os.path.join(self.wdir, self.record_id),
                        str(self.dependency_input["fasta"]),
                    ]
                )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(BuildDatabaseIter.BuildDatabase, BuildDatabaseIter.name, *args, **kwargs)
