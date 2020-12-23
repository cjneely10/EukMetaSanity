from EukMetaSanity import Task, TaskList, program_catch


class RepeatsIter(TaskList):
    """ Task will use NCBI repeats libraries to mask genome.
    Also can generate ab-initio repeat models to mask genome.

    Outputs: mask-fna, mask-tbl, mask-gff3
    Finalizes: mask-tbl

    """
    name = "repeats"
    requires = ["taxonomy"]
    depends = ["repmask.rmout", "repmask.process_repeats"]
    
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-fna": self.input["repmask.rmout"]["mask-fna"],
                "final": [
                    "repmask.process_repeats.rmout",
                    "repmask.process_repeats.rmcat",
                    "repmask.process_repeats.rmtbl",
                    "repmask.rmout.mask-fna",
                    "repmask.rmout.mask-gff3"
                ]
            }
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, RepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
