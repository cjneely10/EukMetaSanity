import os
from pathlib import Path
from typing import List, Union, Type

from yapim import Task, DependencyInput


class Braker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "cds": self.wdir.joinpath(self.record_id + ".cds.fna"),
            "prot": self.wdir.joinpath(self.record_id + ".faa"),
            "gff3": self.wdir.joinpath(self.record_id + ".gff3"),
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run braker
        """
        _tax = []
        if "fungi" in self.input["taxonomy"]["kingdom"]["value"].lower():
            _tax = ["--fungus"]
        _fasta_output = self.input["prots"]
        _added = []
        if not os.path.exists(os.path.join(self.wdir, "GeneMark-ET", "genemark.gtf")) or \
                not os.path.exists(os.path.join(self.wdir, "augustus.hints.gtf")):
            self.parallel(
                self.program[
                    "--cores=%s" % str(self.threads),
                    "--genome=%s" % self.input["fasta"],
                    (*self._get_bams()),
                    "--prot_seq=%s" % _fasta_output,
                    "--workingdir=%s" % self.wdir,
                    (*_tax),
                    "--species=%s" % self.record_id,
                    (*self.added_flags),
                    "--prg=spaln",
                    (*_added)
                ]
            )
        _files = [
            os.path.join(self.wdir, "GeneMark-ET", "genemark.gtf"),
            os.path.join(self.wdir, "augustus.hints.gtf"),
        ]
        _out = os.path.join(self.wdir, "braker.gtf")
        if "--gff3" in self.added_flags:
            _files = [
                os.path.join(self.wdir, "GeneMark-ET", "genemark.gff3"),
                os.path.join(self.wdir, "augustus.hints.gff3"),
            ]
            _out = os.path.join(self.wdir, "braker.gff3")
        # Merge final results if final output failed
        if not os.path.exists(_out):
            self.single(
                self.local["gffread"][
                    (*_files),
                    "--merge", "-G", "-S",
                    "-o", str(self.output["nr_gff3"]),
                    "-g", str(self.input["fasta"]),
                    "-x", str(self.output["cds"]),
                    "-y", str(self.output["prot"]),
                ],
                "30:00"
            )
        # Else just keep
        else:
            self.output["nr_gff3"] = Path(_out)

    def _get_bams(self):
        """ Get list of bam files associated with input as a list in braker-input format

        :return:
        """
        bams = ",".join(list(map(str, self.input["bams"])))
        if len(bams) > 0:
            return ["--bam=%s" % bams]
        return []
