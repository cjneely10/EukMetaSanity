"""
Module holds filtering build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class FilteringIter(TaskList):
    """ TaskList class iterates over filtering tasks

    name: filtering

    requires:

    depends:

    output:

    final:

    """
    name = "filtering"
    requires = ["mapping"]
    depends = [DependencyInput("sambamba.sort", "mapping")]

    class Filtering(Task):
        """
        Task class handles filtering task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "sorted_bams": self.input["sambamba.sort"]["sorted.bams"],
            }

        @program_catch
        def run(self):
            """
            Run filtering
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(FilteringIter.Filtering, FilteringIter.name, *args, **kwargs)
