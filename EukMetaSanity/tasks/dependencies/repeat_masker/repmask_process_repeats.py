import os
from EukMetaSanity import Task, TaskList, program_catch, touch


class ProcessRepeatsIter(TaskList):
    name = "repmask.process_repeats"
    requires = ["taxonomy"]
    depends = ["repmask.repeat_masker"]
    
    class ProcessRepeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "rmout": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.out")
            }
            
        @program_catch
        def run(self):
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            _basename = os.path.basename(str(self.input["root"]["fna"]))
            # Unzip results
            all([
                self.local["gunzip"][os.path.join(rep_dir, "".join((_basename, ".cat.gz")))]()
                for rep_dir in self.input["repmask.repeat_masker"]["libraries"]
                if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat.gz"))))
            ])
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            touch(final_out)
            all([
                (self.local["cat"][os.path.join(rep_dir, "".join((_basename, ".cat")))] >> final_out)()
                for rep_dir in self.input["repmask.repeat_masker"]["libraries"]
                if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat"))))
            ])
            if os.path.getsize(final_out) > 0:
                # Run ProcessRepeats
                self.program[
                    # Input taxonomy from OrthoDB search
                    "-species", self.input["taxonomy"]["taxonomy"].family,
                    "-maskSource", str(self.input["root"]["fna"]),
                    final_out,
                ]()
            
    def __init__(self, *args, **kwargs):
        super().__init__(ProcessRepeatsIter.ProcessRepeats, ProcessRepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
