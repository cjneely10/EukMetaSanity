import os
import glob
import shutil
from typing import List
from pathlib import Path
from EukMetaSanity.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.taxonomy import TaxonomyIter

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            masked_fa_path = os.path.join(self.wdir, "".join((self.record_id, "-mask.out")))
            self.output = [
                masked_fa_path,  # Input FASTA file for ab initio
                masked_fa_path,  # MMseqs database for use in metaeuk
                self.input[0],  # Original input file,
                self.input[2],  # Tax file
            ]

        def run(self) -> None:
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.protocol)(self.input[1])

        # Simple repeat masking using mmseqs
        @program_catch
        def simple(self, input_file: str):
            # Generate the masked sequence
            self.log_and_run(
                self.program_mmseqs[
                    "masksequence",
                    input_file,
                    os.path.join(self.wdir, self.record_id),
                    "--threads", self.threads,
                ]
            )
            # Output as FASTA file
            _fasta_output = "".join((input_file, ".mmseqs_simple.fasta"))
            self.log_and_run(
                self.program_mmseqs[
                    "convert2fasta",
                    os.path.join(self.wdir, self.record_id),
                    _fasta_output,
                ]
            )
            return _fasta_output

        # Complete masking using RepeatModeler/Masker
        def full(self, input_file: str):
            # BuildDatabase and RepeatModeler
            # RepeatMasker and ProcessRepeats
            self._mask(self._model(self.simple(input_file)))

        @program_catch
        def _model(self, input_file: str):
            # Build database
            _name = os.path.join(self.wdir, self.record_id)
            self.log_and_run(
                self.program_builddatabase[
                    "-name", _name,
                    input_file,
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
            return input_file

        @program_catch
        def _mask(self, input_file: str):
            # Perform on de novo results
            de_novo_library = None
            _results_dir = glob.glob("RM_%s*" % str(os.getpid()))
            if len(_results_dir) > 0:
                de_novo_library = os.path.join(_results_dir[0], "consensi.fa")
            # Perform step on each file passed by user
            data_files = []
            if de_novo_library is not None:
                data_files += [de_novo_library]
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
                        input_file,
                    ]
                )
                # Move output file
                if os.path.exists(_search):
                    shutil.move(os.path.dirname(_search), os.path.join(self.wdir))
            # Combine repeat results and process
            self._parse_output(_added_dirs, input_file)

        @program_catch
        def _parse_output(self, repeats_dirs: List[str], input_file: str):
            if len(repeats_dirs) == 0:
                return
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            _basename = os.path.basename(input_file)
            # Unzip results
            all(
                self.log_and_run(self.local["gunzip"][os.path.join(rep_dir, "".join((_basename, ".cat.gz")))])
                for rep_dir in repeats_dirs
            )
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            all(
                self.log_and_run(self.local["cat"][os.path.join(rep_dir, "".join((_basename, ".cat")))] >> final_out)
                for rep_dir in repeats_dirs
            )
            # Run ProcessRepeats
            self.log_and_run(
                self.program_process_repeats[
                    # Input taxonomy from OrthoDB search
                    "-species", TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], float(self.cutoff))[0],
                    "-maskSource", input_file,
                    final_out,
                ]
            )
            # Rename output file
            os.replace(
                input_file + ".masked",
                os.path.join(self.wdir, "".join((self.record_id, "-mask.out")))
            )

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
