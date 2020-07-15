import os
from Bio import SeqIO
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity import Task, TaskList, program_catch


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                os.path.join(self.wdir, self.record_id + ".gff3"),  # Combined results of ab initio + evidence
                os.path.join(self.wdir, self.record_id + ".faa"),  # Proteins
                self.input[4],  # Masked results
            ]

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            # Subset taxonomic database
            subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_db")
            if not os.path.exists(subset_db_outpath):
                self.log_and_run(
                    self.program_mmseqs[
                        "filtertaxseqdb",
                        self.data,
                        subset_db_outpath,
                        "--taxon-list", TaxonomyIter.Taxonomy.get_taxonomy(self.input[3], float(self.cutoff))[1],
                        "--threads", self.threads,
                    ]
                )
            # Run metaeuk
            _outfile = os.path.join(self.wdir, self.record_id)
            self.log_and_run(
                self.program_metaeuk[
                    "easy-predict",
                    self.input[1],
                    subset_db_outpath,
                    _outfile,
                    os.path.join(self.wdir, "tmp"),
                    "--threads", self.threads,
                    "--add-orf-stop",
                ]
            )
            # Convert to GFF3
            self.local["fasta-to-gff3.py"][
                self.input[2], _outfile + ".fas", "-o", os.path.join(self.wdir, "metaeuk.gff3")
            ]()
            # Merge ab initio and initial prediction results
            self.log_and_run(
                self.local["cat"][self.input[0], os.path.join(self.wdir, "metaeuk.gff3")] |
                self.program_gffread[
                    "-o", os.path.join(self.wdir, self.record_id + ".gff3"), "-l", "30",
                    "-y", os.path.join(self.wdir, self.record_id + ".faa"), "-g", self.input[4], "-S",
                    "-Z", "-G", "-M", "-K", "-J", "-Q",
                ]
            )
            # # Rename final output protein sequences
            # out = []
            # for record in SeqIO.parse(_outfile + ".fas", "fasta"):
            #     record.seq = record.seq.upper()
            #     out.append(record)
            # SeqIO.write(out, _outfile + ".faa", "fasta")

    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)


if __name__ == "__main__":
    pass
