"""
Module holds busco build functionality
"""

import os
import glob
from EukMetaSanity import Task, TaskList, prefix
from EukMetaSanity import program_catch, set_complete


class BuscoIter(TaskList):
    """ TaskList class iterates over busco tasks

    name: busco

    requires:

    depends:

    output: results

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
                "results": os.path.join(self.wdir, self.record_id, "short_summary.txt")
            }

        @program_catch
        def run(self):
            """
            Run busco
            """
            results_directory = os.path.dirname(str(self.output["results"]))
            self.parallel(
                self.program[
                    "-i", self.dependency_input["fasta"],
                    "--out_path", results_directory,
                    "-o", prefix(str(self.dependency_input["fasta"])),
                    "-m", self.config["mode"],
                    "-l", self.config["lineage"],
                    "--cpu", self.threads,
                    (*self.added_flags)
                ]
            )
            # Change name of output file
            os.replace(
                glob.glob(os.path.join(self.wdir, self.record_id, "*", "short_summary*.txt"))[0],
                str(self.output["results"])
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(BuscoIter.Busco, BuscoIter.name, *args, **kwargs)
