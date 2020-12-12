import os
from EukMetaSanity import MissingDataError
from EukMetaSanity import Task, TaskList, program_catch


class TaxonomyIter(TaskList):
    """ This class will use `mmseqs` to identify putative taxonomy for an organism.

    Outputs: seq_db, tax_db, tax-report

    Finalizes: tax-report

    """
    name = "taxonomy"
    requires = ["mmseqs.createdb", "mmseqs.taxonomy"]
    
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "seq_db": self.input["mmseqs.createdb"]["db"],
                "tax-report": self.input["mmseqs.taxonomy"]["tax-report"],
                "final": ["tax-report"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
