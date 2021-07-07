import os
from typing import List, Union, Type

from yapim import Task, DependencyInput, touch


class GeneMarkPETAP(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "gtf": self.wdir.joinpath("genemark.gtf"),
            "ab-gff3": self.wdir.joinpath(self.record_id + ".gff3"),
            "prot": self.wdir.joinpath(self.record_id + ".faa")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("GeneMarkProtHint")
        ]

    def run(self):
        """
        Run gmes.petap
        """
        ev_vals = ["--ES"]
        # For now default to 100 intron predictions minimum to use file
        if len(open(str(self.input["GeneMarkProtHint"]["hints"])).readlines()) > 100 and \
                len(open(str(self.input["GeneMarkProtHint"]["evidence"])).readlines()) > 100:
            ev_vals = ["--EP", str(self.input["GeneMarkProtHint"]["hints"]),
                       "--evidence", str(self.input["GeneMarkProtHint"]["evidence"])]
        try:
            if "--ES" not in ev_vals:
                self._run_petap(ev_vals, "gmep.sh")
            else:
                self._run_petap(ev_vals, "gmes.sh")
        except:
            if ev_vals != ["--ES"]:
                self._run_petap(["--ES"], "gmes.sh")
        if os.path.exists(self.output["gtf"]):
            self.single(
                self.local["gffread"][
                    self.output["gtf"], "-G", "-o", str(self.output["ab-gff3"]),
                    "-g", self.input["fasta"],
                    "-y", self.output["prot"]
                ],
                "30:00"
            )
        else:
            touch(str(self.output["ab-gff3"]))
            touch(str(self.output["gtf"]))
            touch(str(self.output["prot"]))

    def _run_petap(self, ev_vals: List[str], script_name: str):
        """ Run gmes_petap.pl

        :param ev_vals: List containing ES- or EP-related command-line flags to pass to run
        """
        script = self.create_script(
            self.program[
                "--sequence", str(self.input["fasta"]),
                (*ev_vals),
                "--cores", self.threads, (*self.added_flags),
                ("--fungus"
                 if self.input["taxonomy"]["kingdom"] is not None and
                    self.input["taxonomy"]["kingdom"]["value"].lower() == "fungi" else "")
            ],
            script_name
        )
        # Run script
        self.parallel(script)
