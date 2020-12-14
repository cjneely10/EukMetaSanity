import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError


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
            tax = self.input["taxonomy"]["taxonomy"].order
            subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % prefix(self.data[0]))
            self.parallel(
                self.program[
                    "filtertaxseqdb",
                    self.data[0],
                    subset_db_outpath,
                    "--taxon-list", tax[1],
                    "--threads", self.threads,
                ],
                "1:00:00"
            )
            # Output as FASTA file
            self.single(
                self.program_mmseqs[
                    "convert2fasta",
                    subset_db_outpath,
                    self.output["fasta"],
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(FilterTaxSeqDBIter.FilterTaxSeqDB, FilterTaxSeqDBIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
