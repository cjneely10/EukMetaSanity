"""
Module holds sambamba.sort build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import program_catch, prefix, set_complete


class SambambaSortIter(TaskList):
    """ TaskList class iterates over sambamba.sort tasks

    name: sambamba.sort

    requires:

    depends: sambamba.view

    output: sorted.bams

    final:

    """
    name = "sambamba.sort"
    requires = []
    depends = [DependencyInput("sambamba.view")]

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
                "sorted.bams": [os.path.join(self.wdir, prefix(db) + ".sorted.bam")
                                for db in self.input["sambamba.view"]["bams"]]
            }

        @program_catch
        def run(self):
            """
            Run sambamba.sort
            """
            for bam_file in self.output["sorted.bams"]:
                out_prefix = os.path.splitext(bam_file)[0]
                self.parallel(
                    self.program[
                        "sort",
                        "-t", self.threads,
                        "-o", out_prefix + ".sorted.bam",
                        bam_file,
                        (*self.added_flags)
                    ]
                )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(SambambaSortIter.SambambaSort, SambambaSortIter.name, *args, **kwargs)
