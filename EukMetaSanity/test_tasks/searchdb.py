import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class MMSeqsSearchDBIter(TaskList):
    class MMSeqsSearchDB(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            
        @program_catch
        def run(self):
            self.local["mmseqs"]["easy-search"][
                self.input["createdb"]["db"]
            ]()
            
    def __init__(self, *args, **kwargs):
        super().__init__(MMSeqsSearchDBIter.MMSeqsSearchDB, "searchdb", ["createdb"],
                         *args, **kwargs)


if __name__ == "__main_":
    pass
