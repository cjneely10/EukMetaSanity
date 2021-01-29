"""
Module holds logic to identify repetitive regions in a genome/MAG
"""
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class RepeatsIter(TaskList):
    """ Task will use NCBI repeats libraries to mask genome.
    Also can generate ab-initio repeat models to mask genome.

    Outputs: mask-fna, mask-tbl, mask-gff3
    Finalizes: mask-tbl

    """
    name = "repeats"
    requires = ["taxonomy"]
    depends = [
        DependencyInput("repmask.rmout"),
        DependencyInput("repmask.process_repeats")
    ]

    class Repeats(Task):
        """
        Run repeats task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "fna": self.input["repmask.rmout"]["mask-fna"],
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
            """
            Run
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Generate task list
        """
        super().__init__(RepeatsIter.Repeats, RepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
