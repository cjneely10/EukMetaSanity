"""
Module holds logic to run `report`-level mmseqs task
"""
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class MMseqsIter(TaskList):
    """ Task searches gene calls through one (or more) mmseqs databases or profile databases

    Outputs: [db-prefix, dynamic]
    Finalizes: [db-prefix, dynamic]

    """
    name = "mmseqs"
    requires = []
    depends = [DependencyInput("mmseqs.convertalis", id_mapping=[("prot", "fasta")])]

    class MMseqs(Task):
        """
        MMseqs task class
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "results": self.input["mmseqs.convertalis"]["results_files"][0],
                "final": ["results"]
            }

        @program_catch
        def run(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(MMseqsIter.MMseqs, MMseqsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
