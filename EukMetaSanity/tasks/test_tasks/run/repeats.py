import os
import glob
import shutil
from typing import List
from pathlib import Path
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity.tasks.test_tasks.run.helpers.taxonomy import get_taxonomy
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class RepeatsIter(TaskList):
    name = "repeats"
    requires = ["taxonomy"]
    
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-fna": os.path.join(self.wdir, self.record_id + "-mask.fna"),
                "mask-tbl": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.tbl"),
                "mask-gff3": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3"),
            }
            
        @program_catch
        def run(self):
            # Mask using NCBI data
            if self.protocol == "simple":
                self.mask(str(self.input["root"]["fna"]))
            # Run ab initio model and incorporate NCBI data
            elif self.protocol == "full":
                self.model()
                self.mask(str(self.input["root"]["fna"]))
            else:
                raise InvalidProtocolError

        @program_catch
        def model(self):
            # Database exists, move on
            if len(glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa.classified"))) > 0:
                return
            _name = os.path.join(self.wdir, self.record_id)
            self.single(
                self.program_builddatabase[
                    "-name", _name,
                    str(self.input["root"]["fna"])
                ]
            )
            script = self.create_script(
                self.program_modeler[
                    "-pa", str(int(self.threads) // 4 or 1),
                    "-engine", "ncbi", "-database", _name
                ],
                "modeler.sh"
            )
            self.parallel(self.local[script])

        @program_catch
        def mask(self, input_file: str):
            data_files = []
            _added_dirs = []
            _file = glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa.classified"))
            if "data" in dir(self):
                data_files += [_file for _file in self.data.split(",") if _file != ""]
            # Perform on optimal taxonomic identification
            data_files += [get_taxonomy(str(self.input["taxonomy"]["tax-report"]), float(self.cutoff), "family")[0]]
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
                # Create contained directory
                self.pm.add_dirs(self.record_id, [_dir])
                _added_dirs.append(self.pm.get_dir(self.record_id, _dir))
                # Call RepeatMasker on modeled repeats in the new directory
                self.parallel(
                    self.program_masker[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", self.pm.get_dir(self.record_id, _dir),
                        input_file,
                    ]
                )
            # Combine repeat results and process
            self.parse_output(_added_dirs, input_file)

        @program_catch
        def parse_output(self, repeats_dirs: List[str], input_file: str):
            if len(repeats_dirs) == 0:
                return
            self.pm.add_dirs(self.record_id, ["repeats_final"])
            _basename = os.path.basename(input_file)
            # Unzip results
            all([
                self.local["gunzip"][os.path.join(rep_dir, "".join((_basename, ".cat.gz")))]()
                for rep_dir in repeats_dirs
                if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat.gz"))))
            ])
            # Combine results into single file
            final_out = os.path.join(self.pm.get_dir(self.record_id, "repeats_final"), "mask.final.cat")
            all([
                (self.local["cat"][os.path.join(rep_dir, "".join((_basename, ".cat")))] >> final_out)()
                for rep_dir in repeats_dirs if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat"))))
            ])
            # Run ProcessRepeats
            self.program_process_repeats[
                # Input taxonomy from OrthoDB search
                "-species", get_taxonomy(str(self.input["taxonomy"]["tax-report"]), float(self.cutoff), "family")[0],
                "-maskSource", input_file,
                final_out,
            ]()
            if os.path.exists(input_file + ".masked"):
                # Rename output file
                os.replace(
                    input_file + ".masked",
                    str(self.output["mask-fna"])
                )
                # Output the repeats file as a gff3 file
                (self.program_rmouttogff3[
                     os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.out")
                 ] > str(self.output["mask-gff3"]))()
            else:
                shutil.copy(input_file, str(self.output["mask-fna"]))
                touch(str(self.output["mask-tbl"]))
                touch(str(self.output["mask-gff3"]))
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, RepeatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
