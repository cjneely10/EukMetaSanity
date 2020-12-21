import os
from EukMetaSanity import Task, TaskList, program_catch, prefix


class ConvertAlisIter(TaskList):
    name = "mmseqs.convertalis"
    requires = []
    depends = ["mmseqs.search", "mmseqs.createdb"]
    
    class ConvertAlis(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results_files": [
                    os.path.join(self.wdir, prefix(db) + ".m8")
                    for db in self.input["mmseqs.search"]["dbs"]]
            }

        @program_catch
        def run(self):
            # Output results
            for data, db, outfile in zip(self.data, self.input["mmseqs.search"]["dbs"], self.output["results_files"]):
                self.parallel(
                    self.program[
                        "convertalis",
                        str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                        data,  # Input augustus-db
                        db,  # Input tax db
                        outfile,  # Output results file
                        "--threads", self.threads,
                        (*self.added_flags)
                    ],
                    "2:00:00"
                )
            
    def __init__(self, *args, **kwargs):
        super().__init__(ConvertAlisIter.ConvertAlis, ConvertAlisIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
