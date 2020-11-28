import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class MMSeqsSearchDBIter(TaskList):
    class MMSeqsSearchDB(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results-file": os.path.join(self.wdir, self.record_id + ".fna"),
                "final": {
                    "results-file": os.path.join(self.wdir, self.record_id + ".fna")
                }
            }
            
        @program_catch
        def run(self):
            out_db = os.path.join(self.wdir, self.record_id + "_db")
            self.log_and_run(
                self.local["mmseqs"]["search"][
                    self.input["createdb"]["db"],
                    self.data,
                    out_db,
                    os.path.join(self.wdir, "tmp"),
                    "--threads", self.threads,
                    (*self.added_flags),
                ]
            )
            self.log_and_run(
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
        super().__init__(MMSeqsSearchDBIter.MMSeqsSearchDB, "searchdb", ["createdb"],
                         *args, **kwargs)


if __name__ == "__main_":
    pass
