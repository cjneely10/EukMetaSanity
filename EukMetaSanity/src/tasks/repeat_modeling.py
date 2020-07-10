import os
import glob
import shutil
from typing import List
from pathlib import Path

from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            masked_db_path = os.path.join(self.wdir, self.record_id + "-mask_db")
            masked_fa_path = masked_db_path[:-3] + ".out"
            self.output = [
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
                self.input[0],  # Original input file,
                self.input[2],  # Tax file
            ]

        def run(self) -> None:
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.protocol)()

        # Simple repeat masking using mmseqs
        @program_catch
        def simple(self):
            # Generate the masked sequence
            self.log_and_run(
                self.program[
                    "masksequence",
                    self.input[1],
                    os.path.join(self.wdir, self.record_id),
                    "--threads", self.threads,
                ]
            )
            # Output as FASTA file
            self.log_and_run(
                self.program[
                    "convert2fasta",
                    os.path.join(self.wdir, self.record_id),
                    self.output[1],
                ]
            )

        # Complete masking using RepeatModeler/Masker
        def full(self):
            # BuildDatabase and RepeatModeler
            self._model()
            # RepeatMasker and ProcessRepeats
            self._mask()

        @program_catch
        def _model(self):
            # Build database
            _name = os.path.join(self.wdir, self.record_id)
            self.log_and_run(
                self.program[
                    "-name", _name,
                    self.input[0],
                ]
            )
            # Run RepeatModeler
            self.log_and_run(
                self.program_modeler[
                    "-pa", self.threads,
                    (*self.added_flags),
                    "-database", _name,
                ]
            )

        @program_catch
        def _mask(self):
            # Perform on de novo results
            de_novo_library = None
            _results_dir = glob.glob(os.path.join(self.wdir, "RM_%s*" % str(os.getpid())))
            if len(_results_dir) > 0:
                de_novo_library = os.path.join(_results_dir[0], "consensi.fa")
            # Perform step on each file passed by user
            data_files = []
            if de_novo_library is not None:
                data_files = [de_novo_library]
            if "data" in dir(self):
                data_files += [_file for _file in self.data.split(",") if _file != ""]
            # Perform on optimal taxonomic identification
            data_files += [TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], float(self.cutoff))[0]]
            _added_dirs = []
            for _search in data_files:
                # Parse for if as file or a RepeatMasker library
                if os.path.exists(str(Path(_search).resolve())):
                    search = ("-lib", _search)
                    _dir = "repeats_" + prefix(_search)
                else:
                    search = ("-species", _search)
                    _dir = "repeats_" + _search.replace(" ", "_")
                # Create contained directory
                self.pm.add_dirs(self.record_id, [_dir])
                _added_dirs.append(self.pm.get_dir(self.record_id, _dir))
                # Call RepeatMasker on modeled repeats in the new directory
                self.log_and_run(
                    self.program_masker[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", self.pm.get_dir(self.record_id, _dir),
                        self.input[0],
                    ]
                )
                # Move output file
                if os.path.exists(_search):
                    shutil.move(_search, os.path.join(self.wdir))
            # Combine repeat results and process
            self._parse_output(_added_dirs)

        @program_catch
        def _parse_output(self, repeats_dirs: List[str]):
            if len(repeats_dirs) == 0:
                return
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            # Unzip results
            all(
                self.log_and_run(self.local["gunzip"]["/".join((rep_dir, "*.cat.gz"))])
                for rep_dir in repeats_dirs
            )
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            self.log_and_run(
                (
                    self.local["cat"][(["/".join((rep_dir, "*.cat")) for rep_dir in repeats_dirs])] >> final_out
                )
            )
            # Run ProcessRepeats
            self.log_and_run(
                self.program_process_repeats[
                    # Input taxonomy from OrthoDB search
                    "-species", TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], float(self.cutoff))[0],
                    "-gff", final_out
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
