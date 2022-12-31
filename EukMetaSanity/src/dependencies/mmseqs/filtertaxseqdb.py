import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, prefix, touch

from EukMetaSanity.mmseqs_taxonomy_report_parser import MMSeqsTaxonomyReportParser


class MMSeqsFilterTaxSeqDB(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        out_results = []
        for i, database in enumerate(self.data):
            if database == "":
                continue
            if "p:" in database:
                database = database[2:]
            db_prefix = prefix(database) + str(i)
            _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
            out_results.append(_outfile)
        self.output = {
            "dbs": out_results,
            "fastas": [out_res + ".fasta" for out_res in out_results],
            "_": self.wdir.joinpath(".done")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run mmseqs.filtertaxseqdb
        """
        assignment = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(self.input["taxonomy"],
                                                                                self.config["level"])
        tax = assignment[1]["taxid"]
        for database, subset_db_outpath, out_fasta in zip(self.data, self.output["dbs"], self.output["fastas"]):
            if tax is None:
                touch(out_fasta)
                continue
            if not os.path.exists(out_fasta):
                self.parallel(
                    self.program[
                        "filtertaxseqdb",
                        database,
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
                    ],
                    "30:00"
                )
        touch(str(self.output["_"]))
