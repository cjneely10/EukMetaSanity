"""
Module holds mmseqs.convertalis build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, DependencyInput, set_complete


class ConvertAlisIter(TaskList):
    """ TaskList class iterates over mmseqs.convertalis tasks

    name: mmseqs.convertalis

    requires:

    depends: mmseqs.search, mmseqs.createdb

    expects:

    output: results_files[List[Path]]

    config:
        mmseqs.convertalis:
          data:
            /path/to/data/odb-mmetsp_db
          program: mmseqs
          FLAGS:
            --format-output query,target,pident,taxid,taxname,taxlineage

    """
    name = "mmseqs.convertalis"
    requires = []
    depends = [
        DependencyInput("mmseqs.search"),
        DependencyInput("mmseqs.createdb")
    ]

    class ConvertAlis(Task):
        """
        Task class handles mmseqs.convertalis task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "results_files": [
                    os.path.join(self.wdir, prefix(database) + "_" + prefix(data) + ".m8")
                    for data, database in zip(self.data, self.input["mmseqs.search"]["dbs"])]
            }

        @program_catch
        def run(self):
            """
            Run mmseqs.convertalis
            """
            # Output results
            for data, database, outfile in zip(self.data, self.input["mmseqs.search"]["dbs"],
                                               self.output["results_files"]):
                if not os.path.exists(outfile):
                    self.parallel(
                        self.program[
                            "convertalis",
                            str(self.input["mmseqs.createdb"]["db"]),  # Input FASTA sequence db
                            data,  # Input augustus-db
                            database,  # Input tax db
                            outfile,  # Output results file
                            "--threads", self.threads,
                            (*self.added_flags)
                        ],
                        "2:00:00"
                    )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(ConvertAlisIter.ConvertAlis, ConvertAlisIter.name, *args, **kwargs)
