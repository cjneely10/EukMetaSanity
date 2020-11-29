import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class MMSeqsSearchDBIter(TaskList):
    name = "searchdb"
    requires = ["createdb"]

    class MMSeqsSearchDB(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results-file": os.path.join(self.wdir, self.record_id + ".m8"),
                "final": {
                    "search": os.path.join(self.wdir, self.record_id + ".m8")
                }
            }
            
        @program_catch
        def run(self):
            out_db = os.path.join(self.wdir, self.record_id + "_db")
            self.run_script(
                self.local["mmseqs"]["search"][
                    self.input["createdb"]["db"],
                    self.data,
                    out_db,
                    os.path.join(self.wdir, "tmp"),
                    "--threads", self.threads,
                    (*self.added_flags),
                ],
                "search.sh"
            )
            self.parallel(
                self.local["mmseqs"]["convertalis"][
                    self.input["createdb"]["db"],
                    self.data,
                    out_db,
                    os.path.join(self.wdir, self.record_id + ".m8"),
                    "--threads", self.threads,
                ],
                "2:00:00"
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(MMSeqsSearchDBIter.MMSeqsSearchDB, MMSeqsSearchDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
