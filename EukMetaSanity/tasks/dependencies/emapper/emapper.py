"""
Module holds emapper build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class EMapperIter(TaskList):
    """ TaskList class iterates over emapper tasks

    name: emapper

    requires:

    depends:

    expects: prot

    output: emapper

    config:
        emapper:
          program: /path/to/emapper.py
          python: /path/to/python2.7
          FLAGS:
            # Provide flags as needed
            -m diamond
    """
    name = "emapper"
    requires = []
    depends = []

    class EMapper(Task):
        """
        Task class handles emapper task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "emapper": os.path.join(self.wdir, self.record_id + ".emapper"),
            }

        @program_catch
        def run(self):
            """
            Run emapper
            """
            self.parallel(
                self.local[self.config["python"]][
                    self.program,
                    "-i", self.dependency_input["prot"],
                    "--output", os.path.join(self.wdir, self.record_id),
                    "--cpu", self.threads,
                    (*self.added_flags),
                ]
            )
            os.replace(
                os.path.join(self.wdir, self.record_id + ".emapper.annotations"),
                str(self.output["emapper"])
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(EMapperIter.EMapper, EMapperIter.name, *args, **kwargs)
