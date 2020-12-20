import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch


class EvidenceIter(TaskList):
    """ Task uses MetaEuk to align protein profiles to genome and output gff file of putative protein locations

    Outputs: metaeuk-gff3, nr-gff3, prot, cds, all_gff3
    Finalizes: metaeuk-gff3, nr-gff3, prot, cds, all_gff3

    """
    name = "evidence"
    requires = ["abinitio.augustus", "abinitio.genemark"]
    depends = ["metaeuk"]
    
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "metaeuk-gff3": os.path.join(self.wdir, "metaeuk.gff3"),  # Metaeuk output
                "nr-gff3": os.path.join(self.wdir, self.record_id + ".nr.gff3"),  # Non-redundant GFF
                "prot": os.path.join(self.wdir, self.record_id + ".faa"),  # NR Proteins
                "cds": os.path.join(self.wdir, self.record_id + ".cds.fna"),  # NR CDS
                "all_gff3": os.path.join(self.wdir, self.record_id + ".all.gff3"),  # Combined gff file
                "final": ["metaeuk-gff3", "nr-gff3", "prot", "cds", "all_gff3"]
            }
            
        @program_catch
        def run(self):
            # Merge final results
            EvidenceIter.Evidence.merge(
                self, [str(self.input["abinitio"]["ab-gff3"]), *out_results],
                str(self.input["root"]["fna"]),
                os.path.join(self.wdir, self.record_id),
            )

        @staticmethod
        def merge(task_object: Task, input_list: List[str], fasta_file: str, out_prefix: str):
            # Convert to gff3 file
            task_object.single(
                task_object.program_gffread[
                    (*input_list), "-G", "--merge",
                    "-o", out_prefix + ".all.gff3"
                ]
            )
            # Replace transcripts with gene identifier and write cds/aa sequences
            task_object.single(
                task_object.local["create-final-annotations.py"][
                    "-f", fasta_file, "-g", out_prefix + ".all.gff3"
                ]
            )
            os.replace(
                out_prefix + ".all.nr.gff3",
                out_prefix + ".nr.gff3",
            )
            os.replace(
                out_prefix + ".all.faa",
                out_prefix + ".faa"
            )
            os.replace(
                out_prefix + ".all.cds.fna",
                out_prefix + ".cds.fna"
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, EvidenceIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
