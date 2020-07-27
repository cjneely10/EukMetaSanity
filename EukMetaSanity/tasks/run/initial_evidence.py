import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            _out = {
                "metaeuk": os.path.join(self.wdir, "metaeuk.gff3"),  # Combined results of ab initio + evidence
                "prot": os.path.join(self.wdir, self.record_id + ".faa"),  # Proteins
                "cds": os.path.join(self.wdir, self.record_id + ".cds.fna"),  # CDS
                "mask": self.input[4],  # Masked results
                "nr_gff3": os.path.join(self.wdir, self.record_id + ".nr.gff3"),  # Non-redundant GFF
                "abinitio": self.input[0],  # Ab initio file
                "tax": self.input[3],  # Taxonomy results file
                "mask_tbl": self.input[5],  # Summarized mask results
            }
            self.output = [
                _out,  # Dictionary for accessing to write final summary
                *list(_out.values()),  # Regular list of values for final path checking
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
            if not os.path.exists(_outfile + ".fas"):
                self.log_and_run(
                    self.program_metaeuk[
                        "easy-predict",
                        self.input[2],
                        subset_db_outpath,
                        _outfile,
                        os.path.join(self.wdir, "tmp"),
                        "--threads", self.threads,
                        "--add-orf-stop",
                    ]
                )
            # Convert to GFF3
            self.log_and_run(
                self.local["fasta-to-gff3.py"][
                    self.input[2], _outfile + ".fas", "-o", os.path.join(self.wdir, "metaeuk.gff3")
                ]
            )
            # Merge to non-redundant set
            self.log_and_run(
                self.local["exonize.py"][
                    "-f", self.input[4],
                    "-g", os.path.join(self.wdir, "metaeuk.gff3"), self.input[0],
                    "-o", os.path.join(self.wdir, self.record_id + ".nr.gff3"),
                    "-c", os.path.join(self.wdir, self.record_id + ".cds.fna"),
                    "-p", os.path.join(self.wdir, self.record_id + ".faa"),
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)


if __name__ == "__main__":
    pass
