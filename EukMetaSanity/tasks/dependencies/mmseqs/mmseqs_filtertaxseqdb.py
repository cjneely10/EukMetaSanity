"""
Module holds mmseqs.filtertaxseqdb build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, set_complete, touch


class FilterTaxSeqDBIter(TaskList):
    """ TaskList class iterates over mmseqs.filtertaxseqdb tasks

    name: mmseqs.filtertaxseqdb

    requires: taxonomy.taxonomy[TaxonomyAssignment]

    depends:

    expects:

    output: dbs[List[Path]], fastas[List[Path]]

    config:
        mmseqs.filtertaxseqdb:
          program: mmseqs
          level: order
          data:
            /path/to/data/odb-mmetsp_db

    """
    name = "mmseqs.filtertaxseqdb"
    requires = ["taxonomy"]
    depends = []

    class FilterTaxSeqDB(Task):
        """
        Task class handles mmseqs.filtertaxseqdb task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
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
                "fastas": [out_res + ".fasta" for out_res in out_results]
            }

        @program_catch
        def run(self):
            """
            Run mmseqs.filtertaxseqdb
            """
            tax = self.input["taxonomy"]["taxonomy"].assignment(self.config["level"], False).tax_id
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

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(FilterTaxSeqDBIter.FilterTaxSeqDB, FilterTaxSeqDBIter.name, *args, **kwargs)
