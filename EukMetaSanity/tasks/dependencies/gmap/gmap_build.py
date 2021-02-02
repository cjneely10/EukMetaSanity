"""
Module holds gmap.build build functionality
"""

import os
from EukMetaSanity import Task, TaskList
from EukMetaSanity import program_catch, set_complete


class GMapBuildIter(TaskList):
    """ TaskList class iterates over gmap.build tasks

    name: gmap.build

    requires:

    depends:

    expects: fasta[Path]

    output: db[Path]

    config:
        gmap.build:
          skip: true
          program: gmapindex

    """
    name = "gmap.build"
    requires = []
    depends = []

    class GMapBuild(Task):
        """
        Task class handles gmap.build task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "db": (os.path.join(self.wdir, self.record_id + "_db") if not self.is_skip else [])
            }

        @program_catch
        def run(self):
            """
            Run gmap.build
            """
            _genome_dir = os.path.dirname(str(self.dependency_input["fasta"]))
            _genome_basename = os.path.basename(str(self.dependency_input["fasta"]))
            self.single(
                self.program[
                    "-d", self.output["db"],
                    "-D", _genome_dir, _genome_basename
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(GMapBuildIter.GMapBuild, GMapBuildIter.name, *args, **kwargs)
