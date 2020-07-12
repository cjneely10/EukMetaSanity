import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.taxonomy import TaxonomyIter


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                os.path.join(self.wdir, self.record_id + ".gff3")  # Combined results of ab initio + evidence
            ]

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            # Subset taxonomic database
            subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_db")
            self.log_and_run(
                self.program_mmseqs[
                    "filtertaxseqdb",
                    self.data,
                    self.input[1],
                    "--taxon-list", TaxonomyIter.Taxonomy.get_taxonomy(self.input[2], float(self.cutoff))[1],
                    subset_db_outpath,
                    "--threads", self.threads,
                ]
            )
            # Run metaeuk
            self.log_and_run(
                self.program_metaeuk[
                    "easy-predict",
                    self.input[1],
                    subset_db_outpath,
                    os.path.join(self.wdir, "tmp"),
                    "--threads", self.threads,
                ]
            )
            # Convert to GFF3
            self.local["fasta-to-gff3.py"](self.input[2], os.path.join(self.wdir, "predsResults.fas"),
                                           "-o", os.path.join(self.wdir, "metaeuk.gff3"))
            # Merge ab initio and initial prediction results
            self.log_and_run(
                self.local["cat"][self.input[0], os.path.join(self.wdir, "metaeuk.gff3")] | self.program_gffread
            )

    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)


if __name__ == "__main__":
    pass
