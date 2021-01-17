"""
Module holds busco build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class BuscoIter(TaskList):
    """ TaskList class iterates over busco tasks

    name: busco

    requires:

    depends:

    output:

    final:

    """
    name = "busco"
    requires = []
    depends = []

    class Busco(Task):
        """
        Task class handles busco task
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
            Run busco
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(BuscoIter.Busco, BuscoIter.name, *args, **kwargs)
