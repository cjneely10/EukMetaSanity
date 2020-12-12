import os
from EukMetaSanity import MissingDataError
from EukMetaSanity import Task, TaskList, program_catch


class TaxonomyIter(TaskList):
    """ This class will use `mmseqs` to identify putative taxonomy for an organism.

    Outputs: seq_db, tax_db, tax-report

    Finalizes: tax-report

    """
    name = "taxonomy"
    requires = []
    
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "seq_db": os.path.join(self.wdir, self.record_id + "_db"),
                "tax-report": os.path.join(self.wdir, "tax-report.txt"),
                "final": ["tax-report"]
            }
            
        @program_catch
        def run(self):
            if not os.path.exists(self.data[0]):
                raise MissingDataError
            # Create mmseqs database
            self.single(
                self.program("mmseqs")[
                    "createdb",
                    self.input["root"]["fna"], self.output["seq_db"]
                ]
            )
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            # Search taxonomy db
            self.parallel(
                self.program("mmseqs")[
                    "taxonomy",
                    self.output["seq_db"],
                    self.data[0],
                    tax_db,
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags("mmseqs")),
                    "--threads", self.threads,
                ]
            )
            # Generate taxonomy report
            self.single(
                self.program("mmseqs")[
                    "taxonomyreport",
                    self.data[0],
                    tax_db,
                    self.output["tax_report"]
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
