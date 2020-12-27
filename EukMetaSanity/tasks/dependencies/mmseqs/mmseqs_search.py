import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, DependencyInput, set_complete


class SearchIter(TaskList):
    name = "mmseqs.search"
    requires = []
    depends = [DependencyInput("mmseqs.createdb")]
    
    class Search(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            outfiles = []
            for db in self.data:
                if db == "":
                    continue
                if "p:" in db:
                    db = db[2:]
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, prefix(db)))
                outfiles.append(_outfile)
            self.output = {
                "dbs": outfiles
            }

        @program_catch
        def run(self):
            for outfile, db_path in zip(self.output["dbs"], self.data):
                if not os.path.exists(outfile + ".index"):
                    # Run search
                    self.parallel(
                        self.program[
                            self.config["subname"],
                            str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                            db_path,  # Input db
                            outfile,  # Output db
                            os.path.join(self.wdir, "tmp"),
                            (*self.added_flags),
                            "--threads", self.threads,
                        ]
                    )
            
    def __init__(self, *args, **kwargs):
        super().__init__(SearchIter.Search, SearchIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
