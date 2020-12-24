from EukMetaSanity import Task, TaskList, program_catch, DependencyInput


class AbInitioAugustusIter(TaskList):
    name = "abinitio.augustus"
    requires = ["repeats"]
    depends = [
        DependencyInput("augustus", "repeats.fna")
    ]
    
    class AbInitioAugustus(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": self.input["augustus"]["ab-gff3"],
                "final": ["augustus.ab-gff3"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioAugustusIter.AbInitioAugustus, AbInitioAugustusIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
