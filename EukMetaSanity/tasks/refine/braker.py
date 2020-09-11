import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter


class BrakerIter(TaskList):
    class Braker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            out = {
                "nr_gff3": None,
            }
            self.output = [
                out,
                *list(out.values()),
            ]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            tax = TaxonomyIter.Taxonomy.get_taxonomy(self.input[3], 0, "kingdom")[0]
            _tax = []
            if "fungi" in tax:
                _tax = ["--fungus"]
            self.log_and_run(
                self.program_braker[
                    "--useexisting",
                    "--cores=%s" % str(self.threads),
                    "--genome=%s" % self.input[0],
                    "--bam=%s" % (",".join((*self.input[-1], *self.input[-2]))),
                    "--prot_seq=%s" % self.input[1],
                    (*_tax),
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(BrakerIter.Braker, "braker", *args, **kwargs)


if __name__ == "__main_":
    pass
