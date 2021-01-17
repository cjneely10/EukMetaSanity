import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, set_complete, DependencyInput


class MetaEukIter(TaskList):
    name = "metaeuk"
    requires = []
    depends = [DependencyInput("mmseqs.filtertaxseqdb")]
    
    class MetaEuk(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "gff3": os.path.join(self.wdir, self.record_id + ".gff3")
            }
            
        @program_catch
        def run(self):
            out_results = []
            for i, db in enumerate(self.input["mmseqs.filtertaxseqdb"]["fastas"]):
                if db == "":
                    continue
                is_profile = []
                if "p:" in db:
                    is_profile.append("--slice-search")
                    db = db[2:]
                db_prefix = prefix(db) + str(i)
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
                if not os.path.exists(_outfile + ".fas"):
                    self.parallel(
                        self.program[
                            "easy-predict",
                            str(self.dependency_input["fasta"]),
                            db,
                            _outfile,
                            os.path.join(self.wdir, "tmp"),
                            "--threads", self.threads,
                            (*self.added_flags),
                            (*is_profile),
                        ]
                    )
                # Convert to GFF3
                self.single(
                    self.local["metaeuk-to-gff3.py"][
                        str(self.dependency_input["fasta"]), _outfile + ".fas", "-o",
                        os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix),
                    ]
                )
                out_results.append(os.path.join(self.wdir, "%s-metaeuk.gff3" % db_prefix))
            self.single(
                self.local["gffread"][
                    (*out_results), "-G", "--cluster-only",
                    "-o", str(self.output["gff3"])
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(MetaEukIter.MetaEuk, MetaEukIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
