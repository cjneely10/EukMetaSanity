import os
from copy import deepcopy
from typing import List, Union, Type

from yapim import Task, DependencyInput


class Braker(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "cds": self.wdir.joinpath(self.record_id + ".cds.fna"),
            "prot": self.wdir.joinpath(self.record_id + ".faa"),
            "gtf": self.wdir.joinpath(self.record_id + ".gtf"),
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
        if "fungi" in self.input["taxonomy"]["kingdom"]["value"].lower():
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
        os.rename(self.wdir.joinpath("augustus.hints.aa"), self.output["prot"])
        os.rename(self.wdir.joinpath("augustus.hints.codingseq"), self.output["cds"])
        os.rename(self.wdir.joinpath("augustus.hints.gtf"), self.output["gtf"])


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
