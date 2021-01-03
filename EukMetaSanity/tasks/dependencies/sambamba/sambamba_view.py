"""
Module holds sambamba.view build functionality
"""

import os
from EukMetaSanity import Task, TaskList
from EukMetaSanity import program_catch, prefix, set_complete


class SambambaViewIter(TaskList):
    """ TaskList class iterates over sambamba.view tasks

    name: sambamba.view

    requires:

    depends:

    expects: sams

    output: bams

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
                "bams": []
            }

        @program_catch
        def run(self):
            """
            Run sambamba.view
            """
            for bam_file in [os.path.join(self.wdir, prefix(db) + ".bam") for db in self.dependency_input["sams"]]:
                self.parallel(
                    self.program[
                        "view",
                        "-S", os.path.splitext(bam_file)[0] + ".sam",
                        "-t", self.threads,
                        "-o", bam_file,
                        (*self.added_flags)
                    ]
                )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(SambambaViewIter.SambambaView, SambambaViewIter.name, *args, **kwargs)
