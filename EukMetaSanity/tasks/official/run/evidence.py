"""
Module holds logic to generate protein-based evidence and to merge together
"""
import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class EvidenceIter(TaskList):
    """ Task uses MetaEuk to align protein profiles to genome and output gff file of putative protein locations

    name: evidence

    requires: abinitio.augustus, abinitio.genemark

    depends: metaeuk.repeats

    output: all_gff3[Path], nr-gff3-tiern[Path], prot-tiern[Path], cds-tiern[Path]

    """
    name = "evidence"
    requires = ["abinitio.augustus", "abinitio.genemark"]
    depends = [DependencyInput("metaeuk", "repeats")]

    class Evidence(Task):
        """
        Evidence class handles merging all lines of evidence (ab initio and protein-based)
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "all_gff3": os.path.join(self.wdir, self.record_id + ".all.gff3"),  # Combined gff file
                "final": ["all_gff3"]
            }
            for i in range(0, 4):
                self.output.update({
                    "nr-gff3-tier%i" % i: os.path.join(self.wdir, self.record_id + ".all.tier%i.nr.gff3" % i),  # NR GFF
                    "prot-tier%i" % i: os.path.join(self.wdir, self.record_id + ".all.tier%i.faa" % i),  # NR Proteins
                    "cds-tier%i" % i: os.path.join(self.wdir, self.record_id + ".all.tier%i.cds.fna" % i),  # NR CDS
                })
                self.output["final"].append("nr-gff3-tier%i" % i)
                self.output["final"].append("prot-tier%i" % i)
                self.output["final"].append("cds-tier%i" % i)

        @program_catch
        def run(self):
            """
            Merge final results
            """
            self.merge(
                [str(self.input["metaeuk"]["gff3"]),
                 self.input["abinitio.augustus"]["ab-gff3"], self.input["abinitio.genemark"]["ab-gff3"]],
                str(self.input["root"]["fasta"]),
                os.path.join(self.wdir, self.record_id),
            )

        def merge(self, input_list: List[str], fasta_file: str, out_prefix: str):
            """  Convert to final non-redundant tiered gff3 files

            :param input_list: List of files to merge
            :param fasta_file: FASTA file with masked data
            :param out_prefix: Output prefix for files
            """
            self.single(
                self.local["gffread"][
                    (*input_list), "-G", "--merge",
                    "-o", out_prefix + ".all.gff3"
                ]
            )
            self.batch([self.local["create-final-annotations.py"][
                            "-f", fasta_file, "-g", out_prefix + ".all.gff3", "-t", i]
                        for i in range(0, 4)])

    def __init__(self, *args, **kwargs):
        """
        Generate task iterator
        """
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)
