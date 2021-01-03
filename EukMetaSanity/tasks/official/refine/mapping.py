"""
Module holds mapping build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class MappingIter(TaskList):
    """ TaskList class iterates over mapping tasks

    name: mapping

    requires:

    depends:

    output:

    final:

    """
    name = "mapping"
    requires = []
    depends = [DependencyInput("gmap"), DependencyInput("hisat2")]

    class Mapping(Task):
        """
        Task class handles mapping task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "sams": [*self.input["gmap"]["sams"], *self.input["hisat2"]["sams"]]
            }

        @program_catch
        def run(self):
            """
            Run mapping
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(MappingIter.Mapping, MappingIter.name, *args, **kwargs)
