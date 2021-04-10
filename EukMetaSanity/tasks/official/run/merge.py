"""
Module holds logic to generate protein-based evidence and to merge together
"""
import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class MergeIter(TaskList):
    """ Task merges all input files from evidence lines into final output

    name: merge

    requires: abinitio.augustus, abinitio.genemark, evidence

    depends:

    output: all_gff3[Path], nr-gff3-tiern[Path], prot-tiern[Path], cds-tiern[Path]

    """
    name = "merge"
    requires = ["abinitio.augustus", "abinitio.genemark", "evidence"]
    depends = []

    class Merge(Task):
        """
        Merge class handles merging all lines of evidence (ab initio and protein-based)
        """

        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "prot": os.path.join(self.wdir, self.record_id + ".faa"),
                "final": ["prot"],
            }

        @program_catch
        def run(self):
            """
            Merge final results
            """
            if os.path.exists(str(self.input["abinitio.genemark"]["ab-gff3"])):
                self.single(
                    self.local["gffread"][
                        "-g", self.input["root"]["fasta"],
                        str(self.input["abinitio.genemark"]["ab-gff3"]),
                        "-y", self.output["prot"]
                    ]
                )

    def __init__(self, *args, **kwargs):
        """
        Generate task iterator
        """
        super().__init__(MergeIter.Merge, MergeIter.name, *args, **kwargs)
