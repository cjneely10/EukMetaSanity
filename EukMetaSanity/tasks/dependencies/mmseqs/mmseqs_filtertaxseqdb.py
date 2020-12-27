import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, set_complete


class FilterTaxSeqDBIter(TaskList):
    name = "mmseqs.filtertaxseqdb"
    requires = ["taxonomy"]
    depends = []
    
    class FilterTaxSeqDB(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            out_results = []
            for db in self.data:
                if db == "":
                    continue
                is_profile = []
                if "p:" in db:
                    is_profile.append("--slice-search")
                    db = db[2:]
                db_prefix = prefix(db)
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
                out_results.append(_outfile)
            self.output = {
                "dbs": out_results,
                "fastas": [out_res + ".fasta" for out_res in out_results]
            }

        @program_catch
        def run(self):
            tax = self.input["taxonomy"]["taxonomy"].assignment(self.config["level"]).tax_id
            for db, subset_db_outpath, out_fasta in zip(self.data, self.output["dbs"], self.output["fastas"]):
                if not os.path.exists(out_fasta):
                    self.parallel(
                        self.program[
                            "filtertaxseqdb",
                            db,
                            subset_db_outpath,
                            "--taxon-list", tax,
                            "--threads", self.threads,
                        ],
                        "1:00:00"
                    )
                    # Output as FASTA file
                    self.single(
                        self.program[
                            "convert2fasta",
                            subset_db_outpath,
                            out_fasta,
                        ]
                    )
            
    def __init__(self, *args, **kwargs):
        super().__init__(FilterTaxSeqDBIter.FilterTaxSeqDB, FilterTaxSeqDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
