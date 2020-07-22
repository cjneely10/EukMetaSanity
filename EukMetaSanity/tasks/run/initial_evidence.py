import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            _out = {
                "gff3": os.path.join(self.wdir, self.record_id + ".gff3"),  # Combined results of ab initio + evidence
                "prot": os.path.join(self.wdir, self.record_id + ".faa"),  # Proteins
                "mask": self.input[4],  # Masked results
                "nr_gff3": os.path.join(self.wdir, self.record_id + ".nr.gff3"),  # Non-redundant GFF
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
            if not os.path.exists(_outfile):
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
            self.local["fasta-to-gff3.py"][
                self.input[2], _outfile + ".fas", "-o", os.path.join(self.wdir, "metaeuk.gff3")
            ]()
            # Merge ab initio and initial prediction results into non-redundant set
            # self.log_and_run(
            #     self.program_bedtools[
            #         "-a", os.path.join(self.wdir, "metaeuk.gff3"),
            #         "-b", self.input[0], "-wa",
            #     ] > os.path.join(self.wdir, self.record_id + ".nr.gff3")
            # )
            # self.log_and_run(
            #     self.local["cat"][self.input[0], os.path.join(self.wdir, "metaeuk.gff3")] |
            #     self.program_gffread[
            #         "-o", os.path.join(self.wdir, self.record_id + ".nr.gff3"), "-S", "-g", self.input[2],
            #         "-Z", "-G", "-M", "-J", "-Q", "-K", "-Y",  # Squash to non-redundant
            #         # "-Z", "-G", "-J", "-M",
            #         "-y", os.path.join(self.wdir, self.record_id + ".tmp.faa")
            #     ]
            # )
            # # Rename proteins
            # record_fp = SeqIO.parse(os.path.join(self.wdir, self.record_id + ".tmp.faa"), "fasta")
            # out = []
            # i = 1
            # for record in record_fp:
            #     record.id = str(record.id) + "_" + str(i)
            #     i += 1
            #     out.append(record)
            # SeqIO.write(out, os.path.join(self.wdir, self.record_id + ".faa"), "fasta")
            # # Remove locus lines
            # self.log_and_run(self.local["sed"]["-i", "/gffcl/d", os.path.join(self.wdir, self.record_id + ".nr.gff3")])
            # Generate complete set, with all redundancies
            self.log_and_run(
                self.local["cat"][self.input[0], os.path.join(self.wdir, "metaeuk.gff3")] |
                self.program_gffread[
                    "-o", os.path.join(self.wdir, self.record_id + ".gff3"), "-g", self.input[2], "-S",
                    "-G", "-M", "--cluster-only", "-J",  # All on top of each other
                ]
            )
            # Create non-redundant gff
            self.log_and_run(
                self.program_gffcompare[
                    os.path.join(self.wdir, self.record_id + ".gff3"), "-o", os.path.join(self.wdir, self.record_id)
                ]
            )
            # Make non-redundant
            self.log_and_run(
                self.program_gffread[
                    os.path.join(self.wdir, self.record_id + ".combined.gtf"), "-G",
                ] > os.path.join(self.wdir, self.record_id + ".nr.gff3")
            )
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/transcript/CDS/g",
                    os.path.join(self.wdir, self.record_id + ".nr.gff3")
                ]
            )
            self.log_and_run(
                self.program_gffread[
                    os.path.join(self.wdir, self.record_id + ".nr.gff3"),
                    "-y", os.path.join(self.wdir, self.record_id + ".faa"),
                    "-g", self.input[2],
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)


if __name__ == "__main__":
    pass
