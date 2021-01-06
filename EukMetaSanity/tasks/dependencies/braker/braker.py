"""
Module holds braker build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class BrakerIter(TaskList):
    """ TaskList class iterates over braker tasks

    name: braker

    requires:

    depends:

    output:

    final:

    """
    name = "braker"
    requires = ["mapping", "taxonomy"]
    depends = [DependencyInput("mmseqs.filtertaxseqdb")]

    class Braker(Task):
        """
        Task class handles braker task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            bams = self._get_bams()
            if len(bams) > 0:
                self.output = {
                    "cds": os.path.join(self.wdir, self.record_id + ".cds.fna"),
                    "prot": os.path.join(self.wdir, self.record_id + ".faa"),
                    "nr_gff3": os.path.join(self.wdir, self.record_id + ".gff3"),
                }
            else:
                self.output = {}

        @program_catch
        def run(self):
            """
            Run braker
            """
            bams = self._get_bams()
            if len(bams) > 0:
                tax = getattr(self.input["taxonomy"]["taxonomy"], self.config["level"])
                _tax = []
                if "fungi" in tax[0]:
                    _tax = ["--fungus"]
                _fasta_output = self.input["mmseqs.filtertaxseqdb"]["fastas"]

        def _get_bams(self):
            """ Get list of bam files associated with input as a list in braker-input format

            :return:
            """
            bams = (",".join(self.input["mapping"]["sorted-bams"]))
            if len(bams) > 0:
                return ["--bam=%s" % ",".join(bams)]
            return []

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(BrakerIter.Braker, BrakerIter.name, *args, **kwargs)
