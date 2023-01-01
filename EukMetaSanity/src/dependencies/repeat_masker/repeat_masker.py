import os
from pathlib import Path
from typing import List, Union, Type

from plumbum import ProcessExecutionError
from yapim import Task, DependencyInput, prefix, touch

from EukMetaSanity.mmseqs_taxonomy_report_parser import MMSeqsTaxonomyReportParser


class RMaskRepeatMasker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        data_files = []
        data_files += [_f for _f in self.data if _f != ""]
        # Perform on optimal taxonomic identification
        assignment = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(self.input["taxonomy"],
                                                                                self.config["level"])
        data_files.append(assignment[1]["value"])
        data_files.append(str(self.input["RModRepeatModeler"]["model"]))
        out = []
        for _search in data_files:
            if "classified" in _search:
                search = ("-lib", str(Path(_search).resolve()))
                _dir = os.path.join(self.wdir, "repeats_" + prefix(_search))
            else:
                search = ("-species", _search)
                _dir = os.path.join(self.wdir, "repeats_" + _search.replace(" ", "_"))
            out.append((search, _dir))
        self.output = {
            "libraries": out,
            "_": self.wdir.joinpath(".done")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("RModRepeatModeler")]

    def run(self):
        """
        Run repmask.repeat_masker
        """
        _added_dirs = []
        for val in self.output["libraries"]:
            search, _dir = val
            # Create contained directory
            if os.path.exists(_dir):
                continue
            os.makedirs(_dir)
            if search[0] == "-lib" and not os.path.exists(search[1]):
                continue
            # Call RepeatMasker on modeled repeats in the new directory
            script = self.create_script(
                self.program[
                    "-pa", self.threads,
                    (*self.added_flags),
                    (*search),
                    "-dir", _dir,
                    self.input["fasta"],
                ],
                "%s.sh" % os.path.basename(_dir)
            )
            try:
                self.parallel(script)
            except ProcessExecutionError:
                continue
        touch(str(self.output["_"]))
