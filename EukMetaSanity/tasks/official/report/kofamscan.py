from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class KoFamScanIter(TaskList):
    """ Task runs kofamscan on a set of gene calls

    Outputs: kegg
    Finalizes: kegg

    """
    name = "kofamscan"
    requires = []
    depends = [DependencyInput("kofamscan.exec_annotation")]
    
    class KoFamScan(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "kegg": self.input["kofamscan.exec_annotation"]["kegg"],
                "final": ["kegg"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(KoFamScanIter.KoFamScan, KoFamScanIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
