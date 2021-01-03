"""
Abstract ab initio task implementing augustus as its search engine
"""

from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class AbInitioAugustusIter(TaskList):
    """
    Class iterates over Augustus Tasks to generate ab initio predicitions

    name: abinitio.augustus

    requires: repeats

    depends: augustus

    output: ab-gff3

    final: ab-gff3
    """
    name = "abinitio.augustus"
    requires = ["repeats"]
    depends = [DependencyInput("augustus", "repeats")]

    class AbInitioAugustus(Task):
        """
        Class calls augustus on input fna
        """
        # @set_complete
        def __init__(self, *args, **kwargs):
            """
            Create Augustus Task
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": self.input["augustus"]["ab-gff3"],
                "final": ["ab-gff3"]
            }

        @program_catch
        def run(self):
            """
            Abstract task is dependency wrapper - no run method is needed
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Create Augustus TaskList
        """
        super().__init__(AbInitioAugustusIter.AbInitioAugustus, AbInitioAugustusIter.name, *args, **kwargs)
