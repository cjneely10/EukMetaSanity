import os
import glob
import shutil
from typing import List
from pathlib import Path
from EukMetaSanity.tasks.utils.helpers import prefix, touch
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            masked_fa_path = os.path.join(self.wdir, "".join((self.record_id, "-mask.fna")))
            masked_db_path = masked_fa_path[:-4] + "_db"
            self.output = [
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
                self.input[0],  # Original input file,
                self.input[2],  # Tax file
                self.input[1],  # MMSeqs db for tax ident
                os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.tbl"),  # Mask results
                os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3")  # Mask results gff3
            ]

        def run_1(self):
            # Call protocol method
            getattr(self, self.protocol)()

        # Simple repeat masking using mmseqs
        @program_catch
        def simple(self):
            # Generate the masked sequence
            self.log_and_run(
                self.program_mmseqs[
                    "masksequence",
                    self.input[1],
                    os.path.join(self.wdir, self.record_id),
                    "--threads", self.threads,
                ]
            )
            _fasta_output = os.path.join(self.wdir, "".join((self.record_id, "-mask_db")))
            # Output as FASTA file
            self.log_and_run(
                self.program_mmseqs[
                    "convert2fasta",
                    os.path.join(self.wdir, self.record_id),
                    _fasta_output,
                ]
            )
            return _fasta_output

        # Complete masking using RepeatModeler/Masker
        def full(self):
            self.simple()
            self._mask(self._model())

        @program_catch
        def _model(self):
            # Build database
            _name = os.path.join(self.wdir, self.record_id)
            self.log_and_run(
                self.program_builddatabase[
                    "-name", _name,
                    self.input[0],
                ]
            )
            new_path = os.path.join(self.wdir, "run.sh")
            self.log_and_run(self.local["cp"][self.local["which"]["run.sh"]().rstrip("\r\n"), self.wdir])
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/g" % (
                        "CMD",
                        " ".join((
                            str(self.program_modeler).replace("/", "\/"),
                            "-pa", str(int(self.threads) // 4 or 1),
                            "-engine", "ncbi", "-database", _name.replace("/", "\/"),
                        ))
                    ), new_path
                ]
            )
            # Run script
            self.log_and_run(self.local[new_path][self.wdir])
            return self.input[0]

        @program_catch
        def _mask(self, input_file: str):
            # Perform on de novo results
            # Perform step on each file passed by user
            data_files = []
            _added_dirs = []
            _file = glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa.classified"))
            if "data" in dir(self):
                data_files += [_file for _file in self.data.split(",") if _file != ""]
            # Perform on optimal taxonomic identification
            data_files += [TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], 0.0, "family")[0]]
            if len(_file) > 0 and os.path.exists(_file[0]):
                data_files.append(_file[0])
            for _search in data_files:
                # Parse for if as file or a RepeatMasker library
                if "RM" in _search:
                    search = ("-lib", str(Path(_search).resolve()))
                    _dir = "repeats_" + prefix(_search)
                else:
                    search = ("-species", _search)
                    _dir = "repeats_" + _search.replace(" ", "_")
                # Do not repeat if step is already present
                if os.path.exists(os.path.join(os.path.dirname(self.wdir), _dir)):
                    continue
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
            # Combine repeat results and process
            self._parse_output(_added_dirs, input_file)

        @program_catch
        def _parse_output(self, repeats_dirs: List[str], input_file: str):
            if len(repeats_dirs) == 0:
                return
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            _basename = os.path.basename(input_file)
            # Unzip results
            all([
                self.log_and_run(self.local["gunzip"][os.path.join(rep_dir, "".join((_basename, ".cat.gz")))])
                for rep_dir in repeats_dirs
                if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat.gz"))))
            ])
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            all([
                self.log_and_run(self.local["cat"][os.path.join(rep_dir, "".join((_basename, ".cat")))] >> final_out)
                for rep_dir in repeats_dirs
            ])
            # Run ProcessRepeats
            self.log_and_run(
                self.program_process_repeats[
                    # Input taxonomy from OrthoDB search
                    "-species", TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], 0.0, "family")[0],
                    "-maskSource", input_file,
                    final_out,
                ]
            )
            if os.path.exists(input_file + ".masked"):
                # Rename output file
                os.replace(
                    input_file + ".masked",
                    os.path.join(self.wdir, "".join((self.record_id, "-mask.fna")))
                )
                # Output the repeats file as a gff3 file
                self.log_and_run(
                    self.program_rmouttogff3[
                        os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.out")
                    ] > os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3")
                )
            else:
                touch(os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3"))

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
