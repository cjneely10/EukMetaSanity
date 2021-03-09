"""
Module holds logic to identify repetitive regions in a genome/MAG
"""
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class RepeatsIter(TaskList):
    """ TaskList class iterates over repeats tasks

    name: repeats

    requires: taxonomy

    depends: repmask.rmout, repmask.process_repeats

    output: fna[Path]

    final: repmask.{process_repeats.rmtbl,rmout.mask-fna,rmout.mask-gff3}[Path]

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
                "mask-fna": self.input["repmask.rmout"]["mask-fna"],
                "final": [
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
