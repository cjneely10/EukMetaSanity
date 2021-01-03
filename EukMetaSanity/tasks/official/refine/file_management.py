"""
Module holds file_management build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class FileManagementIter(TaskList):
    """ TaskList class iterates over file_management tasks

    name: file_management

    requires:

    depends:

    output:

    final:

    """
    name = "file_management"
    requires = ["transcriptomes", "rnaseq"]
    depends = [DependencyInput("sambamba.sort", "transcriptomes"),
               DependencyInput("sambamba.sort", "rnaseq")]

    class FileManagement(Task):
        """
        Task class handles file_management task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "bams": [*self.input["sambamba.sort"]["transcriptomes"], *self.input["sambamba.sort"]["rnaseq"]]
            }

        @program_catch
        def run(self):
            """
            Run file_management
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(FileManagementIter.FileManagement, FileManagementIter.name, *args, **kwargs)
