"""
Module holds transcriptomes build functionality
"""

from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import program_catch, set_complete


class TranscriptomeIter(TaskList):
    """ TaskList class iterates over transcriptomes tasks

    name: transcriptomes

    requires:

    depends: gmap

    output: sams

    final:

    """
    name = "transcriptomes"
    requires = []
    depends = [DependencyInput("gmap")]

    class Transcriptome(Task):
        """
        Task class handles transcriptomes task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "sams": self.input["gmap"]["sams"]
            }

        @program_catch
        def run(self):
            """
            Run transcriptomes
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(TranscriptomeIter.Transcriptome, TranscriptomeIter.name, *args, **kwargs)
