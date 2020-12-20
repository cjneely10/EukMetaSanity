from EukMetaSanity import Task, TaskList, program_catch


class KoFamScanIter(TaskList):
    """ Task runs kofamscan on a set of gene calls

    Outputs: kegg
    Finalizes: kegg

    """
    name = "kofamscan"
    requires = []
    depends = ["kofamscan.exec_annotation"]
    
    class KoFamScan(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "final": ["kofamscan.exec_annotation.kegg"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, KoFamScanIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
