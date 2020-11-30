import os
from EukMetaSanity import Task, TaskList, program_catch


class KoFamScanIter(TaskList):
    """ Task runs kofamscan on a set of gene calls

    Outputs: kegg
    Finalizes: kegg

    """
    name = "kofamscan"
    requires = []
    
    class KoFamScan(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "kegg": os.path.join(self.wdir, self.record_id + ".kegg"),
                "final": ["kegg"]
            }
            
        @program_catch
        def run(self):
            self.parallel(
                self.program[
                    "--cpu", self.threads,
                    "--format", "detail",
                    (*self.added_flags),
                    "-o", self.output["kegg"],
                    "--tmp-dir", os.path.join(self.wdir, "tmp"),
                    self.input["root"]["prot"],
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, KoFamScanIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
