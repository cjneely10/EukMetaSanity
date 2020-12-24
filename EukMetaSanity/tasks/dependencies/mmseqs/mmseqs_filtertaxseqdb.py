import os
from EukMetaSanity import Task, TaskList, program_catch, prefix


class FilterTaxSeqDBIter(TaskList):
    name = "mmseqs.filtertaxseqdb"
    requires = ["taxonomy"]
    depends = []
    
    class FilterTaxSeqDB(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "fasta": os.path.join(self.wdir, self.record_id + ".faa")
            }

        @program_catch
        def run(self):
            tax = self.input["taxonomy"]["taxonomy"].order.tax_id
            subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % prefix(self.data[0]))
            self.parallel(
                self.program[
                    "filtertaxseqdb",
                    self.data[0],
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
