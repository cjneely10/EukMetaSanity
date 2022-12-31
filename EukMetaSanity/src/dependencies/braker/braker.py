import os
from copy import deepcopy
from typing import List, Union, Type

from yapim import Task, DependencyInput, touch

from EukMetaSanity.mmseqs_taxonomy_report_parser import MMSeqsTaxonomyReportParser


class Braker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "possible_files": {
                "cds": self.wdir.joinpath(self.record_id + ".cds.fna"),
                "prot": self.wdir.joinpath(self.record_id + ".faa"),
                "gtf": self.wdir.joinpath(self.record_id + ".gtf"),
                "joined": self.wdir.joinpath("braker.gtf")
            },
            "_": self.wdir.joinpath(".done")
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
        tax = []
        assignment = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(self.input["taxonomy"], "kingdom")
        if "fungi" in assignment[1]["value"].lower():
            tax = ["--fungus"]
        self.parallel(
            self.program[
                "--useexisting",
                "--cores=%s" % str(self.threads),
                "--genome=%s" % self.input["fasta"],
                (*self._get_bams()),
                "--prot_seq=%s" % self.input["prots"],
                "--workingdir=%s" % self.wdir,
                (*tax,),
                "--species=%s" % self.record_id,
                (*self._filter_provided_flags()),
            ]
        )
        for key, file in [("prot", self.wdir.joinpath("augustus.hints.aa")),
                          ("cds", self.wdir.joinpath("augustus.hints.codingseq")),
                          ("gtf", self.wdir.joinpath("augustus.hints.gtf"))]:
            if file.exists():
                os.rename(file, self.output["possible_files"][key])
        touch(str(self.output["_"]))

    def _filter_provided_flags(self) -> [str]:
        """ Remove flags that would affect program behaviour

        :return: List of valid user-provided flags
        """
        flags = deepcopy(self.added_flags)
        arg = "--skipGetAnnoFromFasta"
        if arg in flags:
            flags.remove(arg)
        arg = "--gff3"
        if arg in flags:
            flags.remove(arg)
        return flags

    def _get_bams(self):
        """ Get list of bam files associated with input as a list in braker-input format

        :return:
        """
        bams = ",".join(list(map(str, self.input["bams"])))
        if len(bams) > 0:
            return ["--bam=%s" % bams]
        return []
