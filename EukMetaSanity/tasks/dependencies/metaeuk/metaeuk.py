"""
Module holds metaeuk build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, set_complete, DependencyInput


class MetaEukIter(TaskList):
    """ TaskList class iterates over metaeuk tasks

    name: metaeuk

    requires:

    depends: mmseqs.filtertaxseqdb

    expects: fasta

    output: gff3

    config:
        metaeuk:
          program: metaeuk
          data:
            /path/to/data/odb-mmetsp_db
          # Pass any flags to metaeuk required
          FLAGS:
            --min-length 30
            --metaeuk-eval 0.0001
            --split-memory-limit 12G
            -s 7
            --cov-mode 0
            -c 0.3
            -e 100
            --max-overlap 0

    """
    name = "metaeuk"
    requires = []
    depends = [DependencyInput("mmseqs.filtertaxseqdb")]

    class MetaEuk(Task):
        """
        Task class handles metaeuk task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "gff3": os.path.join(self.wdir, self.record_id + ".evidence.gff3")
            }

        @program_catch
        def run(self):
            """
            Run metaeuk
            """
            out_results = []
            for i, database in enumerate(self.input["mmseqs.filtertaxseqdb"]["fastas"]):
                if database == "":
                    continue
                is_profile = []
                if "p:" in database:
                    is_profile.append("--slice-search")
                    database = database[2:]
                db_prefix = prefix(database) + str(i)
                _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
                if not os.path.exists(_outfile + ".fas"):
                    self.parallel(
                        self.program[
                            "easy-predict",
                            str(self.dependency_input["fasta"]),
                            database,
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
        """
        Instantiate TaskList
        """
        super().__init__(MetaEukIter.MetaEuk, MetaEukIter.name, *args, **kwargs)
