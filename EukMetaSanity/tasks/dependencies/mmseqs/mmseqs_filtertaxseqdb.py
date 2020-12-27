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
            # TODO: Remove hard-coding of order-level search
            tax = self.input["taxonomy"]["taxonomy"].order.tax_id
            for db in self.data:
                subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % prefix(db))
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
                        self.output["fasta"],
                    ]
                )
                for file in os.listdir(self.wdir):
                    file = os.path.join(self.wdir, file)
                    if subset_db_outpath in file:
                        self.local["rm"][file]()
            
    def __init__(self, *args, **kwargs):
        super().__init__(FilterTaxSeqDBIter.FilterTaxSeqDB, FilterTaxSeqDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
