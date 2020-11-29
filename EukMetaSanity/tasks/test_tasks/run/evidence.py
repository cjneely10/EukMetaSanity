import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity.tasks.test_tasks.run.helpers.taxonomy import get_taxonomy
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


class EvidenceIter(TaskList):
    name = "evidence"
    requires = ["repeats", "taxonomy", "abinitio"]
    
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
                    self.parallel(
                        self.program_mmseqs[
                            "filtertaxseqdb",
                            db,
                            subset_db_outpath,
                            "--taxon-list", get_taxonomy(
                                str(self.input["root"]["fna"]), 0, self.level,  # Allow for level override by user
                            )[1],
                            "--threads", self.threads,
                        ],
                        "30:00"
                    )
                # Run metaeuk
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
                if not os.path.exists(_outfile + ".fas"):
                    self.parallel(
                        self.program_metaeuk[
                            "easy-predict",
                            str(self.input["repeats"]["mask-fna"]),
                            subset_db_outpath,
                            _outfile,
                            os.path.join(self.wdir, "tmp"),
                            "--threads", self.threads,
                            (*self.added_flags),
                            (*is_profile),
                        ]
                    )
                # Convert to GFF3
                self.single(
                    self.local["fasta-to-gff3.py"][
                        str(self.input["root"]["fna"]), _outfile + ".fas", "-o",
                        os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix),
                    ]
                )
                out_results.append(os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix))
            self.single(
                self.program_gffread[
                    (*out_results), "-G", "--cluster-only",
                    "-o", str(self.output["metaeuk-gff3"])
                ]
            )
            # Merge final results
            EvidenceIter.Evidence.merge(
                self, [str(self.input["abinitio"]["ab-gff3"]), *out_results],
                str(self.input["root"]["fna"]),
                os.path.join(self.wdir, self.record_id),
            )

        @staticmethod
        def merge(task_object: Task, input_list: List[str], fasta_file: str, out_prefix: str):
            # Convert to gff3 file
            task_object.program_gffread[
                (*input_list), "-G", "--merge",
                "-o", out_prefix + ".all.gff3"
            ]()
            # Replace transcripts with gene identifier and write cds/aa sequences
            task_object.local["create-final-annotations.py"][
                "-f", fasta_file, "-g", out_prefix + ".all.gff3"
            ]()
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
