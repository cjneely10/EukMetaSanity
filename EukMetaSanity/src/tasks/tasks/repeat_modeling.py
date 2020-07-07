import os
from typing import List
from pathlib import Path
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.helpers import prefix
from EukMetaSanity.bin.fastagff3_to_gb import write_genbank
from EukMetaSanity.src.tasks.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.utils.config_manager import ConfigManager
from EukMetaSanity.src.tasks.task_class import Task, TaskList, program_catch

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            masked_db_path = os.path.join(self.wdir, self.record_id + "-mask_db")
            masked_fa_path = masked_db_path[:-3] + ".out"
            self.output = {Data.Type.OUT: [
                self.input[Data.Type.IN][0],  # Taxonomic results for ab initio prediction
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
            ]}
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.config[ConfigManager.PROTOCOL])()

        # Simple repeat masking using mmseqs
        @program_catch
        def simple(self):
            masked_db_path = self.output[Data.Type.OUT][2]
            # Generate the masked sequence
            self.log_and_run(
                self.program[
                    "masksequence",
                    self.input[Data.Type.IN][2],
                    masked_db_path,
                    "--threads", self.threads,
                ]
            )
            # Output as FASTA file
            self.log_and_run(
                self.program[
                    "convert2fasta",
                    masked_db_path,
                    self.output[Data.Type.OUT][1],
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
            self.log_and_run(
                self.program[
                    "-name", self.output[Data.Type.OUT][2],
                    (*self.added_flags),
                    self.input[Data.Type.IN][1],
                ]
            )
            # Run RepeatModeler
            self.log_and_run(
                self.program2[
                    "-pa", self.threads,
                    (*self.added_flags),
                    "-database", self.output[Data.Type.OUT][2],
                ]
            )

        @program_catch
        def _mask(self):
            # Perform step on each file
            data_files = [_file for _file in self.config[ConfigManager.DATA].split(",") if _file != ""]
            _added_dirs = []
            for _search in data_files:
                # Parse for if as file or a RepeatMasker library
                if os.path.exists(str(Path(_search).resolve())):
                    search = ("-lib", _search)
                else:
                    search = ("-species", _search)
                _dir = "repeats_" + prefix(_search)
                # Create contained directory
                self.pm.add_dirs(self.record_id, [_dir])
                _added_dirs.append(self.pm.get_dir(self.record_id, _dir))
                # Call RepeatMasker on modeled repeats in the new directory
                self.log_and_run(
                    self.program3[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", self.pm.get_dir(self.record_id, _dir),
                        "-database", self.output[Data.Type.OUT][2],
                    ]
                )
            # Combine repeat results and process
            self._parse_output(_added_dirs)

        @program_catch
        def _parse_output(self, repeats_dirs: List[str]):
            # Unzip results
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            all(
                self.log_and_run(self.local["gunzip", "/".join((rep_dir, "*.cat.gz"))])
                for rep_dir in repeats_dirs
            )
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            self.log_and_run(
                (
                    self.local[
                        "cat", (*["/".join((rep_dir, "*.cat")) for rep_dir in repeats_dirs])
                    ] >> final_out
                )
            )
            # Run ProcessRepeats
            self.log_and_run(
                self.program4[
                    "-species", TaxonomyIter.Taxonomy.get_taxonomy(self.input[Data.Type.IN][0]), final_out
                ]
            )
            # Create GFF3
            out_gff3 = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.gff3")
            self.log_and_run(
                self.program5[
                    os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.out"),
                ] > out_gff3
            )
            # Write as genbank
            write_genbank(
                self.input[Data.Type.IN][1],
                out_gff3,
                self.output[Data.Type.OUT][1]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
