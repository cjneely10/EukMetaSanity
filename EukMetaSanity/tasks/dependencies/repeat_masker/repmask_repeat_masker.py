import os
import shutil
from typing import List
from EukMetaSanity import ProcessExecutionError
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch


class RepeatMaskerIter(TaskList):
    name = "repmask.repeat_masker"
    requires = ["taxonomy", "repmod.repeat_modeler"]
    
    class RepeatMasker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "mask-fna": os.path.join(self.wdir, self.record_id + "-mask.fna"),
                "mask-tbl": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.tbl"),
                "mask-gff3": os.path.join(os.path.dirname(self.wdir), "repeats_final", "mask.final.gff3"),
            }
            
        @program_catch
        def run(self):
            pass

        @program_catch
        def mask(self, input_file: str):
            data_files = []
            _added_dirs = []
            _file = str(self.input["repmod.repeat_modeler"]["model"])
            if "data" in dir(self):
                data_files += [_file for _file in self.data if _file != ""]
            # Perform on optimal taxonomic identification
            data_files += [get_taxonomy(str(self.input["taxonomy"]["tax-report"]), 5.0, "family")[0]]
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
                script = self.create_script(
                    self.program_masker[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", self.pm.get_dir(self.record_id, _dir),
                        input_file,
                    ],
                    "mask.sh"
                )
                try:
                    self.parallel(script)
                except ProcessExecutionError as e:
                    continue
            # Combine repeat results and process
            self.parse_output(_added_dirs, input_file)

        @program_catch
        def parse_output(self, repeats_dirs: List[str], input_file: str):
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
            touch(final_out)
            all([
                (self.local["cat"][os.path.join(rep_dir, "".join((_basename, ".cat")))] >> final_out)()
                for rep_dir in repeats_dirs if os.path.exists(os.path.join(rep_dir, "".join((_basename, ".cat"))))
            ])
            if os.path.getsize(final_out) > 0:
                # Run ProcessRepeats
                self.program_process_repeats[
                    # Input taxonomy from OrthoDB search
                    "-species", get_taxonomy(str(self.input["taxonomy"]["tax-report"]), float(self.cutoff), "family")[
                        0],
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
        super().__init__(RepeatMaskerIter.RepeatMasker, RepeatMaskerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
