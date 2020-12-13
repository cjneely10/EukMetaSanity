import os
import shutil
from EukMetaSanity import Task, TaskList, program_catch, touch


class RepeatMaskerOutIter(TaskList):
    name = "repmask.rmout"
    requires = ["repmask.process_repeats"]
    
    class RepeatModelerOut(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-fna": os.path.join(self.wdir, self.record_id + "-mask.fna"),
                "mask-tbl": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.tbl"),
                "mask-gff3": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3"),
            }
            
        @program_catch
        def run(self):
            if os.path.exists(str(self.input["root"]["fna"]) + ".masked"):
                # Rename output file
                os.replace(
                    str(self.input["root"]["fna"]) + ".masked",
                    str(self.output["mask-fna"])
                )
                # Output the repeats file as a gff3 file
                (self.program[
                     os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.out")
                 ] > str(self.output["mask-gff3"]))()
            else:
                shutil.copy(str(self.input["root"]["fna"]), str(self.output["mask-fna"]))
                touch(str(self.output["mask-tbl"]))
                touch(str(self.output["mask-gff3"]))
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerOutIter.RepeatModelerOut, RepeatMaskerOutIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
