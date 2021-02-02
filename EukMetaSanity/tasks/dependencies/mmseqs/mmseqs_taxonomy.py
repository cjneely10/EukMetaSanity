"""
Module holds mmseqs.taxonomy build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class TaxonomyIter(TaskList):
    """ TaskList class iterates over mmseqs.taxonomy tasks

    name: mmseqs.taxonomy

    requires:

    depends: mmseqs.createdb

    output: tax-report[Path]

    config:
        mmseqs.taxonomy:
          program: mmseqs
          data:
            /path/to/data/odb-mmetsp_db
          # Pass any flags to mmseqs required
          FLAGS:
            --remove-tmp-files
            -s 7
            --min-seq-id 0.40
            -c 0.3
            --cov-mode 0
            --split-memory-limit 12G

    """
    name = "mmseqs.taxonomy"
    requires = []
    depends = [DependencyInput("mmseqs.createdb")]

    class Taxonomy(Task):
        """
        Task class handles mmseqs.taxonomy task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "tax-report": os.path.join(self.wdir, "tax-report.txt"),
            }

        @program_catch
        def run(self):
            """
            Run mmseqs.taxonomy
            """
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            # Search taxonomy db
            self.parallel(
                self.program[
                    "taxonomy",
                    self.input["mmseqs.createdb"]["db"],
                    self.data[0],
                    tax_db,
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags),
                    "--threads", self.threads,
                ]
            )
            # Generate taxonomy report
            self.single(
                self.program[
                    "taxonomyreport",
                    self.data[0],
                    tax_db,
                    self.output["tax-report"]
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)
