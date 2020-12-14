import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class ConvertAlisIter(TaskList):
    name = "mmseqs.convertalis"
    requires = []
    depends = ["mmseqs.linsearch", "mmseqs.createdb"]
    
    class ConvertAlis(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results_file": os.path.join(self.wdir, self.record_id + ".m8")
            }
            
        @program_catch
        def run(self):
            # Output results
            self.parallel(
                self.program[
                    "convertalis",
                    str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                    self.data[0],  # Input augustus-db
                    self.input["mmseqs.linsearch"]["db"],  # Input tax db
                    self.output["results-file"],  # Output results file
                    "--threads", self.threads
                ],
                "2:00:00"
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(ConvertAlisIter.ConvertAlis, ConvertAlisIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
