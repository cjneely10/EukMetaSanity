"""
Module holds sambamba.view build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class SambambaViewIter(TaskList):
    """ TaskList class iterates over sambamba.view tasks

    name: sambamba.view

    requires:

    depends:

    output:

    final:

    """
    name = "sambamba.view"
    requires = []
    depends = []

    class SambambaView(Task):
        """
        Task class handles sambamba.view task
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
            Run sambamba.view
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(SambambaViewIter.SambambaView, SambambaViewIter.name, *args, **kwargs)
