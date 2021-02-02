"""
Module holds combined-evidence build functionality
"""

from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import program_catch, set_complete


class EvidenceIter(TaskList):
    """ TaskList class iterates over combined-evidence tasks

    name: combined-evidence

    requires: mapping, filtering, taxonomy

    depends: braker

    output: cds[Path], prot[Path], nr_gff3[Path]

    final: cds[Path], prot[Path], nr_gff3[Path]

    """
    name = "combined-evidence"
    requires = ["mapping", "filtering", "taxonomy"]
    depends = [DependencyInput("braker")]

    class Evidence(Task):
        """
        Task class handles combined-evidence task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "cds": self.input["braker"]["cds"],
                "prot": self.input["braker"]["prot"],
                "nr_gff3": self.input["braker"]["nr_gff3"],
                "final": ["cds", "prot", "nr_gff3"]
            }

        @program_catch
        def run(self):
            """
            Run combined-evidence
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)
