import os
import shutil

from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput, set_complete


class RepeatMaskerOutIter(TaskList):
    name = "repmask.rmout"
    requires = []
    depends = [DependencyInput("repmask.process_repeats")]

    class RepeatModelerOut(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-gff3": os.path.join(self.wdir, "mask.final.gff3"),
                "fna": os.path.join(self.wdir, self.record_id + ".mask.fna")
            }

        @program_catch
        def run(self):
            input_file = str(self.dependency_input) + ".masked"
            if os.path.exists(input_file):
                os.replace(
                    input_file,
                    str(self.output["mask-fna"])
                )
                # Output the repeats file as a gff3 file
                self.single(
                    (self.program[
                         self.input["repmask.process_repeats"]["rmout"]
                     ] > str(self.output["mask-gff3"]))
                )
            else:
                touch(str(self.output["mask-gff3"]))
                shutil.copyfile(
                    str(self.dependency_input),
                    str(self.output["mask-fna"])
                )

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerOutIter.RepeatModelerOut, RepeatMaskerOutIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
