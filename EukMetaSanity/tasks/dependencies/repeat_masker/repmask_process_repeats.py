import os
from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput


class ProcessRepeatsIter(TaskList):
    name = "repmask.process_repeats"
    requires = ["taxonomy"]
    depends = [DependencyInput("repmask.repeat_masker")]

    class ProcessRepeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "rmout": os.path.join(self.wdir, "mask.final.out"),
                "rmtbl": os.path.join(self.wdir, "mask.final.tbl"),
                "rmcat": os.path.join(self.wdir, "mask.final.cat"),
            }

        @program_catch
        def run(self):
            _basename = os.path.basename(str(self.input["root"]["fna"]))
            cat_files = []
            for rep_dir in self.input["repmask.repeat_masker"]["libraries"]:
                _file = os.path.join(rep_dir[1], "".join((_basename, ".cat.gz")))
                if os.path.exists(_file) and os.path.getsize(_file) > 0:
                    cat_files.append(_file)
            # Unzip results
            all([
                self.local["gunzip"][_file]()
                for _file in cat_files if os.path.exists(_file)
            ])
            cat_files = []
            for rep_dir in self.input["repmask.repeat_masker"]["libraries"]:
                _file = os.path.join(rep_dir[1], "".join((_basename, ".cat")))
                if os.path.exists(_file) and os.path.getsize(_file) > 0:
                    cat_files.append(_file)
            # Combine results into single file
            final_out = os.path.join(self.wdir, "mask.final.cat")
            touch(final_out)
            all([
                (self.local["cat"][os.path.splitext(_file)[0]] >> final_out)()
                for _file in cat_files
                if os.path.exists(os.path.splitext(_file)[0]) and os.path.getsize(os.path.splitext(_file)[0]) > 0
            ])
            if os.path.getsize(final_out) > 0:
                # Run ProcessRepeats
                self.single(
                    self.program[
                        # Input taxonomy from OrthoDB search
                        "-species", self.input["taxonomy"]["taxonomy"].family.value,
                        "-maskSource", str(self.input["root"]["fna"]),
                        final_out,
                    ]
                )
            else:
                touch(str(self.output["rmcat"]))
                touch(str(self.output["rmtbl"]))
                touch(str(self.output["rmout"]))

    def __init__(self, *args, **kwargs):
        super().__init__(ProcessRepeatsIter.ProcessRepeats, ProcessRepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
