import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, DependencyInput, set_complete


class ConvertAlisIter(TaskList):
    name = "mmseqs.convertalis"
    requires = []
    depends = [
        DependencyInput("mmseqs.search"),
        DependencyInput("mmseqs.createdb")
    ]
    
    class ConvertAlis(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results_files": [
                    os.path.join(self.wdir, prefix(db) + "_" + prefix(data) + ".m8")
                    for data, db in zip(self.data, self.input["mmseqs.search"]["dbs"])]
            }

        @program_catch
        def run(self):
            # Output results
            for data, db, outfile in zip(self.data, self.input["mmseqs.search"]["dbs"], self.output["results_files"]):
                if not os.path.exists(outfile):
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
