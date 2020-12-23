import os
from EukMetaSanity import Task, TaskList, program_catch, touch


class RepeatMaskerOutIter(TaskList):
    name = "repmask.rmout"
    requires = []
    depends = ["repmask.process_repeats"]

    class RepeatModelerOut(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-gff3": os.path.join(self.wdir, "mask.final.gff3"),
            }

        @program_catch
        def run(self):
            if os.path.getsize(str(self.input["repmask.process_repeats"]["rmout"])) > 0:
                # Output the repeats file as a gff3 file
                (self.program[
                     self.input["repmask.process_repeats"]["rmout"]
                 ] > str(self.output["mask-gff3"]))()
            else:
                touch(str(self.output["mask-gff3"]))

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerOutIter.RepeatModelerOut, RepeatMaskerOutIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
