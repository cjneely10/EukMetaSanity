"""

"""
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class AbInitioGeneMarkIter(TaskList):
    name = "abinitio.genemark"
    requires = ["taxonomy", "repeats"]
    depends = [DependencyInput("gmes.gffread", "repeats")]
    
    class AbInitioGeneMark(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": self.input["gmes.gffread"]["ab-gff3"],
                "final": ["gmes.gffread.ab-gff3"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioGeneMarkIter.AbInitioGeneMark, AbInitioGeneMarkIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
