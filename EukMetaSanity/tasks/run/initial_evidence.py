import os
from typing import List
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter

"""
Add protein evidence using MetaEuk
Provides custom summary for easy in use in next pipelines

"""


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            _merged_out = os.path.join(self.wdir, self.record_id + ".nr.gff3")
            if os.path.exists(_merged_out):
                os.remove(_merged_out)
            _out = {
                "metaeuk": os.path.join(self.wdir, "metaeuk.gff3"),  # Combined results of ab initio + evidence
                "prot": os.path.join(self.wdir, self.record_id + ".faa"),  # Proteins
                "cds": os.path.join(self.wdir, self.record_id + ".cds.fna"),  # CDS
                "mask": self.input[4],  # Masked results
                "nr_gff3": os.path.join(self.wdir, self.record_id + ".nr.gff3"),  # Non-redundant GFF
                "abinitio": self.input[0],  # Ab initio file
                "tax": self.input[3],  # Taxonomy results file
                "mask_tbl": self.input[5],  # Summarized mask results
                "mask_gff3": self.input[6],  # Mask gff3 file
                "fna": self.input[2],  # Original fna file
                "all_gff": os.path.join(self.wdir, self.record_id + ".all.gff3"),  # Combined gff file
            }
            self.output = [
                _out,  # Dictionary for accessing to write final summary
                *list(_out.values()),  # Regular list of values for final path checking
            ]

        @program_catch
        def run_1(self):
            # Subset taxonomic database
            out_results = []
            for db in self.data.split(","):
                if db == "":
                    continue
                is_profile = []
                if "p:" in db:
                    is_profile.append("--slice-search")
                    db = db[2:]
                db_prefix = prefix(db)
                subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % db_prefix)
                if not os.path.exists(subset_db_outpath):
                    self.log_and_run(
                        self.program_mmseqs[
                            "filtertaxseqdb",
                            db,
                            subset_db_outpath,
                            "--taxon-list", TaxonomyIter.Taxonomy.get_taxonomy(
                                self.input[3], 0, self.level,  # Allow for level override by user
                            )[1],
                            "--threads", self.threads,
                        ]
                    )
                # Run metaeuk
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
                if not os.path.exists(_outfile + ".fas"):
                    self.log_and_run(
                        self.program_metaeuk[
                            "easy-predict",
                            self.input[4],
                            subset_db_outpath,
                            _outfile,
                            os.path.join(self.wdir, "tmp"),
                            "--threads", self.threads,
                            (*self.added_flags),
                            (*is_profile),
                        ]
                    )
                # Convert to GFF3
                self.log_and_run(
                    self.local["fasta-to-gff3.py"][
                        self.input[2], _outfile + ".fas", "-o", os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix),
                    ]
                )
                out_results.append(os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix))
            self.log_and_run(
                self.program_gffread[
                    (*out_results), "-G", "--cluster-only",
                    "-o", os.path.join(self.wdir, "metaeuk.gff3")
                ]
            )
            # Merge final results
            EvidenceIter.Evidence.merge(
                self, [self.input[0], *out_results],
                self.input[2],
                os.path.join(self.wdir, self.record_id),
            )
            
        @staticmethod
        def merge(task_object: Task, input_list: List[str], fasta_file: str, out_prefix: str):
            # Convert to gff3 file
            task_object.log_and_run(
                task_object.program_gffread[
                    (*input_list), "-G", "--merge",
                    "-o", out_prefix + ".all.gff3"
                ]
            )
            task_object.local["sed"][
                "-i", '/gffcl/d',
                out_prefix + ".all.gff3"
            ]()
            # Replace transcripts with gene identifier and write cds/aa sequences
            task_object.log_and_run(
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
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)


if __name__ == "__main__":
    pass
