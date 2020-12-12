import os
from EukMetaSanity import Task, TaskList, program_catch


class CreateDBIter(TaskList):
    """ This class will create an MMseqs database

    """
    name = "mmseqs.createdb"
    requires = []

    class MMSeqsCreateDB(Task):
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
                    self.input["root"]["fna"],
                    self.output["db"]
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(CreateDBIter.MMSeqsCreateDB, CreateDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
