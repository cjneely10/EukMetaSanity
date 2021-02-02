"""
Module holds hisat2 build functionality
"""
import os
from pathlib import Path
from typing import List, Tuple
from EukMetaSanity import Task, TaskList, MissingDataError, DependencyInput
from EukMetaSanity import program_catch, prefix, set_complete


class Hisat2Iter(TaskList):
    """ TaskList class iterates over Hisat2 tasks

    name: hisat2

    requires:

    depends: hisat2.build

    expects:

    output: sams[Path]

    config:
        hisat2:
          skip: true
          program: hisat2
          # Should be in format (excluding spaces around tab):
          # file-basename \t /path/to/r1.fq[,/path/to/r2.fq][;/path/to/r3.fq[,/path/to/r4.fq]]
          rnaseq: /path/to/rnaseq-mapping-file
          FLAGS:
            ""

    """
    name = "hisat2"
    requires = []
    depends = [DependencyInput("hisat2.build")]

    class Hisat2(Task):
        """
        Task class handles Hisat2 task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            if self.is_skip:
                sams = []
            else:
                sams = [os.path.join(self.wdir, prefix(pair[0])) for pair in self.get_rna_read_pairs()]
            self.output = {
                "sams": sams
            }

        @program_catch
        def run(self):
            """
            Run hisat2
            """
            rna_pairs = self.get_rna_read_pairs()
            for pair in rna_pairs:
                out_prefix = os.path.join(self.wdir, prefix(pair[0]))
                if len(pair) > 1:
                    _parse_args = ["-1", pair[0], "-2", pair[1]]
                else:
                    _parse_args = ["-U", pair[0]]
                # Align
                self.parallel(
                    self.program[
                        "-p", self.threads,
                        "-x", self.input["hisat2.build"]["db"],
                        (*_parse_args),
                        "-S", out_prefix + ".sam",
                        (*self.added_flags),
                    ]
                )

        def get_rna_read_pairs(self) -> List[Tuple[str, ...]]:
            """ Parse rna_seq file into list of rna pairs to analyze

            :return: List of read pairs
            """
            if self.developer_mode:
                return []
            _path = str(Path(self.config["rnaseq"]).resolve())
            if not os.path.exists(_path):
                raise MissingDataError("Input rna_seq mapping file not found!")
            file_ptr = open(_path, "r")
            for line in file_ptr:
                if line.startswith(self.record_id):
                    pairs_string = line.rstrip("\r\n").split("\t")[1].split(";")
                    return [tuple(pair.split(",")) for pair in pairs_string]
            return []

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(Hisat2Iter.Hisat2, Hisat2Iter.name, *args, **kwargs)
