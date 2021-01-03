"""
Module holds sambamba.sort build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class SambambaSortIter(TaskList):
    """ TaskList class iterates over sambamba.sort tasks

    name: sambamba.sort

    requires:

    depends:

    output:

    final:

    """
    name = "sambamba.sort"
    requires = []
    depends = []

    class SambambaSort(Task):
        """
        Task class handles sambamba.sort task
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
            Run sambamba.sort
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(SambambaSortIter.SambambaSort, SambambaSortIter.name, *args, **kwargs)
