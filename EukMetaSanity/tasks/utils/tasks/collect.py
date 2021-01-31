"""
Module holds collect build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


# TODO: Implement
class CollectInputIter(TaskList):
    """ TaskList class iterates over collect tasks

    name: collect

    requires:

    depends:

    output:

    final:

    """
    name = "collect"
    requires = []
    depends = []

    class CollectInput(Task):
        """
        Task class handles collect task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                
            }

        @program_catch
        def run(self):
            """
            Run collect
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(CollectInputIter.CollectInput, CollectInputIter.name, *args, **kwargs)
