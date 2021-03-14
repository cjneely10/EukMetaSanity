"""
Module holds repmask.rmout build functionality
"""
import os
import shutil
from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput, set_complete


class RepeatMaskerOutIter(TaskList):
    """ TaskList class iterates over repmask.rmout tasks

    name: repmask.rmout

    requires:

    depends: repmask.process_repeats

    expects: fasta[Path]

    output: mask-gff3[Path], mask-fna[Path]

    config:
        repmask.rmout:
          program: rmOutToGFF3.pl

    """
    name = "repmask.rmout"
    requires = []
    depends = [DependencyInput("repmask.process_repeats")]

    class RepeatModelerOut(Task):
        """
        Task class handles repmask.rmout task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-gff3": os.path.join(self.wdir, "mask.final.gff3"),
                "mask-fna": os.path.join(self.wdir, self.record_id + ".mask.fna")
            }

        @program_catch
        def run(self):
            """
            Run repmask.rmout
            """
            input_file = str(self.dependency_input["fasta"]) + ".masked"
            if os.path.exists(input_file):
                os.replace(
                    input_file,
                    str(self.output["mask-fna"])
                )
            else:
                shutil.copyfile(
                    str(self.dependency_input["fasta"]),
                    str(self.output["mask-fna"])
                )
            # Output the repeats file as a gff3 file
            self.single(
                (self.program[
                     self.input["repmask.process_repeats"]["rmout"]
                 ] > str(self.output["mask-gff3"])),
                "3:00:00"
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(RepeatMaskerOutIter.RepeatModelerOut, RepeatMaskerOutIter.name, *args, **kwargs)
