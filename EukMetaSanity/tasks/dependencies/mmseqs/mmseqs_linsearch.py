import os
from EukMetaSanity import Task, TaskList, program_catch


class LinSearchIter(TaskList):
    name = "mmseqs.linsearch"
    requires = []
    depends = ["mmseqs.createdb"]
    
    class LinSearch(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results_file": os.path.join(self.wdir, self.record_id + ".m8")
            }
            
        @program_catch
        def run(self):
            # Run search
            out_db = os.path.join(self.wdir, self.record_id + "_db")
            self.parallel(
                self.program[
                    "linsearch",
                    str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                    self.data[0],  # Input augustus-db
                    out_db,  # Output tax db
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags),
                    "--threads", self.threads,
                ]
            )
            # Output results
            self.parallel(
                self.program_mmseqs[
                    "convertalis",
                    str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                    self.data[0],  # Input augustus-db
                    out_db,  # Input tax db
                    self.output["results-file"],  # Output results file
                    "--threads", self.threads,
                    "--format-output", "query,target,pident,taxid,taxname,taxlineage",
                ],
                "2:00:00"
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(LinSearchIter.LinSearch, LinSearchIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
