import os
from EukMetaSanity import MissingDataError
from EukMetaSanity import Task, TaskList, program_catch


class TaxonomyIter(TaskList):
    name = "mmseqs.taxonomy"
    requires = ["mmseqs.createdb"]
    
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "tax-report": os.path.join(self.wdir, "tax-report.txt"),
            }
            
        @program_catch
        def run(self):
            if not os.path.exists(self.data[0]):
                raise MissingDataError
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
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
