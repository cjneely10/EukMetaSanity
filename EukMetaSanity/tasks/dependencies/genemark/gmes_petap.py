"""
Module holds gmes.petap build functionality
"""
import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete, touch


class GeneMarkPetapIter(TaskList):
    """ TaskList class iterates over gmes.petap tasks

    name: gmes.petap

    requires: taxonomy.taxonomy[TaxonomyAssignment]

    depends: gmes.prothint

    expects: fasta[Path]

    output: gtf[Path]

    config:
        gmes.petap:
          program: gmes_petap.pl
          FLAGS:
            --min_contig 100
            --max_contig 100000000
            --max_gap 5000
            --max_mask 5000
            --min_contig_in_predict 100
            --min_gene_in_predict 10
            --gc_donor 0.001
            --max_intron 10000
            --max_intergenic 50000
            --soft_mask auto

    """
    name = "gmes.petap"
    requires = ["taxonomy"]
    depends = [DependencyInput("gmes.prothint")]

    class GeneMarkPetap(Task):
        """
        Task class handles gmes.petap task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "gtf": os.path.join(self.wdir, "genemark.gtf"),
                "ab-gff3": os.path.join(self.wdir, self.record_id + ".gmes.gff3")
            }

        @program_catch
        def run(self):
            """
            Run gmes.petap
            """
            ev_vals = ["--ES"]
            # For now default to 100 intron predictions minimum to use file
            if len(open(str(self.input["gmes.prothint"]["hints"])).readlines()) > 100 and \
                    len(open(str(self.input["gmes.prothint"]["evidence"])).readlines()) > 100:
                ev_vals = ["--EP", str(self.input["gmes.prothint"]["hints"]),
                           "--evidence", str(self.input["gmes.prothint"]["evidence"])]
            try:
                self._run_petap(ev_vals)
            except:
                if ev_vals != ["--ES"]:
                    self._run_petap(["--ES"])
            if not os.path.exists(str(self.output["ab-gff3"])):
                self._run_petap(["--ES"])
            if os.path.exists(self.output["gtf"]):
                self.single(
                    self.local["gffread"][
                        self.output["gtf"], "-G", "-o", str(self.output["ab-gff3"])
                    ],
                    "30:00"
                )
            else:
                touch(str(self.output["ab-gff3"]))
                touch(self.output["gtf"])

        def _run_petap(self, ev_vals: List[str]):
            """ Run gmes_petap.pl

            :param ev_vals: List containing ES- or EP-related command-line flags to pass to run
            """
            script = self.create_script(
                self.program[
                    "--sequence", str(self.dependency_input["fasta"]),
                    (*ev_vals),
                    "--cores", self.threads, (*self.added_flags),
                    ("--fungus"
                     if self.input["taxonomy"]["taxonomy"].kingdom is not None and
                        self.input["taxonomy"]["taxonomy"].kingdom.value.lower() == "fungi" else "")
                ],
                "abinitio.sh"
            )
            # Run script
            self.parallel(script)

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(GeneMarkPetapIter.GeneMarkPetap, GeneMarkPetapIter.name, *args, **kwargs)
