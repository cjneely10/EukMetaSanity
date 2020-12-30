import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class CreateDBIter(TaskList):
    """ This class will create an MMseqs database

    name: mmseqs.createdb

    requires: None

    output keys: db

    finalizes: None

    """
    name = "mmseqs.createdb"
    requires = []
    depends = []

    class MMSeqsCreateDB(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "db": os.path.join(self.wdir, self.record_id + "_db")
            }

        @program_catch
        def run(self):
            self.single(
                self.program[
                    "createdb",
                    self.dependency_input["fasta"],
                    self.output["db"]
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(CreateDBIter.MMSeqsCreateDB, CreateDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
