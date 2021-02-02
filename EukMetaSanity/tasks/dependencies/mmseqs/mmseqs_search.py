"""
Module holds mmseqs.search build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, DependencyInput, set_complete


class SearchIter(TaskList):
    """ TaskList class iterates over mmseqs.search tasks

    name: mmseqs.search

    requires:

    depends: mmseqs.createdb

    expects:

    output: dbs[List[Path]]

    config:
        mmseqs.search:
          data:
            /path/to/data/odb-mmetsp_db
          program: mmseqs
          subname: search  # Can use `linsearch` if linear index is present for database
          FLAGS:
            --split-memory-limit 12G
            --cov-mode 0
            -c 0.6
            -e 0.01
            --remove-tmp-files

    """
    name = "mmseqs.search"
    requires = []
    depends = [DependencyInput("mmseqs.createdb")]

    class Search(Task):
        """
        Task class handles mmseqs.search task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            outfiles = []
            for i, database in enumerate(self.data):
                if database == "":
                    continue
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, prefix(database) + str(i)))
                outfiles.append(_outfile)
            self.output = {
                "dbs": outfiles
            }

        @program_catch
        def run(self):
            """
            Run mmseqs.search
            """
            for outfile, db_path in zip(self.output["dbs"], self.data):
                if not os.path.exists(outfile + ".index"):
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
        """
        Instantiate TaskList
        """
        super().__init__(SearchIter.Search, SearchIter.name, *args, **kwargs)
