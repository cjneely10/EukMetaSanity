from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class MMseqsIter(TaskList):
    """ Task searches gene calls through one (or more) mmseqs databases or profile databases

    Outputs: [db-prefix, dynamic]
    Finalizes: [db-prefix, dynamic]

    """
    name = "mmseqs"
    requires = []
    depends = [DependencyInput("mmseqs.search")]
    
    class MMseqs(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "final": ["mmseqs.search.dbs"]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(MMseqsIter.MMseqs, MMseqsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
