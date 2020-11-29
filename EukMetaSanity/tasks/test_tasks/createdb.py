import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class MMSeqsCreateDBIter(TaskList):
    name = "createdb"
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
                self.local["mmseqs"]["createdb"][self.input["root"]["fna"], self.output["db"]]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(MMSeqsCreateDBIter.MMSeqsCreateDB, MMSeqsCreateDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
