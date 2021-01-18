"""
Module holds quality build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class QualityIter(TaskList):
    """ TaskList class iterates over quality tasks

    name: quality

    requires:

    depends:

    output:

    final:

    """
    name = "quality"
    requires = []
    depends = []

    class Quality(Task):
        """
        Task class handles quality task
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
            Run quality
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(QualityIter.Quality, QualityIter.name, *args, **kwargs)
