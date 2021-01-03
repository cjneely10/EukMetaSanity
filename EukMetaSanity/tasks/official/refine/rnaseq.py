"""
Module holds rnaseq build functionality
"""

from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import program_catch, set_complete


class RNASeqIter(TaskList):
    """ TaskList class iterates over rnaseq tasks

    name: rnaseq

    requires:

    depends: hisat2

    output:

    final:

    """
    name = "rnaseq"
    requires = []
    depends = [DependencyInput("hisat2")]

    class RNASeq(Task):
        """
        Task class handles rnaseq task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "sams": self.input["hisat2"]["sams"]
            }

        @program_catch
        def run(self):
            """
            Run rnaseq
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(RNASeqIter.RNASeq, RNASeqIter.name, *args, **kwargs)
