import os
from EukMetaSanity import Task, TaskList, program_catch, prefix


class MetaEukIter(TaskList):
    name = "metaeuk"
    requires = ["repeats"]
    depends = ["mmseqs.filtertaxseqdb"]
    
    class MetaEuk(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "gff3": os.path.join(self.wdir, self.record_id + ".gff3")
            }
            
        @program_catch
        def run(self):
            db = self.data[0]
            db_prefix = prefix(db)
            _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
            self.parallel(
                self.program[
                    "easy-predict",
                    str(self.input["repeats"]["mask-fna"]),
                    str(self.input["mmseqs.filtertaxseqdb"]["db"]),
                    _outfile,
                    os.path.join(self.wdir, "tmp"),
                    "--threads", self.threads,
                    (*self.added_flags),
                ]
            )
            # Convert to GFF3
            self.single(
                self.local["metaeuk-to-gff3.py"][
                    str(self.input["root"]["fna"]), _outfile + ".fas", "-o",
                    os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix),
                ]
            )
            self.single(
                self.local["gffread"][
                    os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix), "-G", "--cluster-only",
                    "-o", str(self.output["gff3"])
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(MetaEukIter.MetaEuk, MetaEukIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
