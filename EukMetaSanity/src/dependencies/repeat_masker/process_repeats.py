import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, touch


class RMaskProcessRepeats(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)
        self.output = {
            "rmout": os.path.join(self.wdir, "mask.final.out"),
            "rmtbl": os.path.join(self.wdir, "mask.final.tbl"),
            "rmcat": os.path.join(self.wdir, "mask.final.cat"),
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("RMaskRepeatMasker")
        ]

    def run(self):
        """
        Run repmask.process_repeats
        """
        _basename = os.path.basename(str(self.input["fasta"]))
        cat_files = []
        for rep_dir in self.input["RMaskRepeatMasker"]["libraries"]:
            _file = os.path.join(rep_dir[1], "".join((_basename, ".cat.gz")))
            if os.path.exists(_file) and os.path.getsize(_file) > 0:
                cat_files.append(_file)
        # Unzip results
        all([
            self.local["gunzip"][_file]()
            for _file in cat_files if os.path.exists(_file)
        ])
        cat_files = []
        for rep_dir in self.input["RMaskRepeatMasker"]["libraries"]:
            _file = os.path.join(rep_dir[1], "".join((_basename, ".cat")))
            if os.path.exists(_file) and os.path.getsize(_file) > 0:
                cat_files.append(_file)
        # Combine results into single file
        final_out = os.path.join(self.wdir, "mask.final.cat")
        if os.path.exists(final_out):
            os.remove(final_out)
        touch(final_out)
        all([(self.local["cat"][_file] >> final_out)() for _file in cat_files])
        if os.path.getsize(final_out) > 0:
            # Run ProcessRepeats
            self.single(
                self.program[
                    # Input taxonomy from OrthoDB search
                    "-species", self.input["taxonomy"]["family"]["value"],
                    "-maskSource", str(self.input["fasta"]),
                    (*self.added_flags),
                    final_out,
                ]
            )
        touch(str(self.output["rmcat"]))
        touch(str(self.output["rmtbl"]))
        touch(str(self.output["rmout"]))
