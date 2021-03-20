"""
Module holds evidence build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class EvidenceIter(TaskList):
    """ TaskList class iterates over evidence tasks

    name: evidence

    requires:

    depends:

    expects:

    output:

    final:

    config:

    """
    name = "evidence"
    requires = []
    depends = [DependencyInput("metaeuk")]

    class Evidence(Task):
        """
        Task class handles evidence task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "prot-gff3": self.input["metaeuk"]["gff3"],
                "prot": self.input["metaeuk"]["prot"]
            }

        @program_catch
        def run(self):
            """
            Run evidence
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)
