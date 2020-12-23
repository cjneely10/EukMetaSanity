from EukMetaSanity import Task, TaskList, program_catch


class AbInitioGeneMarkIter(TaskList):
    name = "abinitio.genemark"
    requires = ["taxonomy", "repeats"]
    depends = ["gmes.gffread"]
    
    class AbInitioGeneMark(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "final": ["gmes.gffread.ab-gff3"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioGeneMarkIter.AbInitioGeneMark, AbInitioGeneMarkIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
