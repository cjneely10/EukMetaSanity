import os
from EukMetaSanity import Task, TaskList, program_catch, data_catch


class SearchIter(TaskList):
    name = "mmseqs.linsearch"
    requires = []
    depends = ["mmseqs.createdb"]
    
    class LinSearch(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "db": os.path.join(self.wdir, self.record_id + "_db")
            }

        @data_catch
        @program_catch
        def run(self):
            # Run search
            self.parallel(
                self.program[
                    self.config["subname"],
                    str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                    self.data[0],  # Input db
                    self.output["results_file"],  # Output db
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags),
                    "--threads", self.threads,
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(SearchIter.LinSearch, SearchIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
