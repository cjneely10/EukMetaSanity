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

    expects: fasta[Path]

    output: gff3[Path], prot[Path]

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
    depends = []

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
                "gff3": os.path.join(self.wdir, self.record_id + ".evidence.gff3"),
                "prot": os.path.join(self.wdir, self.record_id + ".evidence.faa")
            }

        # TODO: Revert to multiple dbs
        @program_catch
        def run(self):
            """
            Run metaeuk
            """
            if len(self.data) == 0:
                return
            database = self.data[0]
            is_profile = []
            if "p:" in database:
                is_profile.append("--slice-search")
                database = database[2:]
            db_prefix = prefix(database)
            _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
            # Run MetaEuk
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
            # Write in GFF3 format
            self.single(
                self.local["metaeuk-to-gff3.py"][
                    str(self.dependency_input["fasta"]), _outfile + ".fas", "-o",
                    str(self.output["gff3"]),
                ],
                time_override="8:00",
                memory_override="6G"
            )
            # Rename output file
            os.replace(_outfile + ".fas", str(self.output["prot"]))

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(MetaEukIter.MetaEuk, MetaEukIter.name, *args, **kwargs)
