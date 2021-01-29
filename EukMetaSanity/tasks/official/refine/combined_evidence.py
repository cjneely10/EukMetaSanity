"""
Module holds combined-evidence build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class EvidenceIter(TaskList):
    """ TaskList class iterates over combined-evidence tasks

    name: combined-evidence

    requires:

    depends:

    output:

    final:

    """
    name = "combined-evidence"
    requires = []
    depends = []

    class Evidence(Task):
        """
        Task class handles combined-evidence task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {}

        @program_catch
        def run(self):
            """
            Run combined-evidence
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)
