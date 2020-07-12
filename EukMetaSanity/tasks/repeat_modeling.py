import os
import shutil
import datetime
from typing import List
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
            masked_db_path = masked_fa_path[:-4] + "_db"
            self.output = [
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
                self.input[0],  # Original input file,
                self.input[2],  # Tax file
                self.input[1],  # MMSeqs db for tax ident
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
        def full(self, input_file: str):
            # BuildDatabase and RepeatModeler
            # RepeatMasker and ProcessRepeats
            self._mask(*self._model(self.simple(input_file)))

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
            _now = RepeatsIter.Repeats.roundTime(datetime.datetime.now())
            # _now = datetime.datetime.now()
            print(_now)
            # Run RepeatModeler
            self.log_and_run(
                self.program_modeler[
                    "-pa", self.threads,
                    (*self.added_flags),
                    "-database", _name,
                ]
            )
            return input_file, _now

        @program_catch
        def _mask(self, input_file: str, _recorded_start_time: datetime.datetime):
            # Perform on de novo results
            from dask.distributed import get_worker
            # Perform step on each file passed by user
            data_files = []
            _file = RepeatsIter.Repeats._get_results_file(_recorded_start_time)
            if _file is not None:
                data_files.append(_file)
            if "data" in dir(self):
                data_files += [_file for _file in self.data.split(",") if _file != ""]
            # Perform on optimal taxonomic identification
            data_files += [TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], float(self.cutoff))[0]]
            _added_dirs = []
            for _search in data_files:
                # Parse for if as file or a RepeatMasker library
                if os.path.exists(_search):
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

        @staticmethod
        def _get_results_file(_time: datetime.datetime):
            # Get list of files to search
            _formatted_time = _time.strftime("%A%b%d%H%M%S%Y")
            print(_formatted_time)
            _files = [_file for _file in os.listdir(os.getcwd()) if "RM" in _file]
            if len(_files) == 1:
                return os.path.join(_files[0], "consensi.fa.classified")
            else:
                _possible_time = (_time + datetime.timedelta(0, 1)).strftime("%A%b%d%H%M%S%Y")
                print(_possible_time)
                for _file in _files:
                    if _formatted_time in _file or _possible_time in _file:
                        return os.path.join(_file, "consensi.fa.classified")

        @staticmethod
        def roundTime(dt=None, roundTo=1):
            """Round a datetime object to any time lapse in seconds
            dt : datetime.datetime object, default now.
            roundTo : Closest number of seconds to round to, default 1 minute.
            Author: Thierry Husson 2012 - Use it as you want but don't blame me.
            """
            if dt == None: dt = datetime.datetime.now()
            seconds = (dt.replace(tzinfo=None) - dt.min).seconds
            rounding = (seconds + roundTo / 2) // roundTo * roundTo
            return dt + datetime.timedelta(0, rounding - seconds, -dt.microsecond)

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
