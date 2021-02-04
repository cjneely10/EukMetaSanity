"""
Module holds repmask.process_repeats build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, touch, DependencyInput, set_complete


class ProcessRepeatsIter(TaskList):
    """ TaskList class iterates over repmask.process_repeats tasks

    name: repmask.process_repeats

    requires: taxonomy.taxonomy[TaxonomyAssignment]

    depends: repmask.repeat_masker

    expects: fasta[Path]

    output: rmout[Path], rmtbl[Path], rmcat[Path]

    config:
        repmask.process_repeats:
          program: ProcessRepeats
          FLAGS:
            -nolow

    """
    name = "repmask.process_repeats"
    requires = ["taxonomy"]
    depends = [DependencyInput("repmask.repeat_masker")]

    class ProcessRepeats(Task):
        """
        Task class handles repmask.process_repeats task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "rmout": os.path.join(self.wdir, "mask.final.out"),
                "rmtbl": os.path.join(self.wdir, "mask.final.tbl"),
                "rmcat": os.path.join(self.wdir, "mask.final.cat"),
            }

        @program_catch
        def run(self):
            """
            Run repmask.process_repeats
            """
            _basename = os.path.basename(str(self.dependency_input["fasta"]))
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
            if os.path.exists(final_out):
                os.remove(final_out)
            touch(final_out)
            all([(self.local["cat"][_file] >> final_out)() for _file in cat_files])
            if os.path.getsize(final_out) > 0 and os.path.exists(str(self.output["rmtbl"])):
                # Run ProcessRepeats
                self.single(
                    self.program[
                        # Input taxonomy from OrthoDB search
                        "-species", self.input["taxonomy"]["taxonomy"].family.value,
                        "-maskSource", str(self.dependency_input["fasta"]),
                        (*self.added_flags),
                        final_out,
                    ]
                )
            else:
                touch(str(self.output["rmcat"]))
                touch(str(self.output["rmtbl"]))
                touch(str(self.output["rmout"]))

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(ProcessRepeatsIter.ProcessRepeats, ProcessRepeatsIter.name, *args, **kwargs)
