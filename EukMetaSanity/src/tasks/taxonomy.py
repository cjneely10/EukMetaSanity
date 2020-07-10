import os
from typing import Tuple
from EukMetaSanity import Task, TaskList, program_catch

"""
Determine the taxonomy of the Eukaryotic MAG

"""


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            seq_db = os.path.join(self.wdir, self.record_id + "_db")
            # Expected output
            self.output = [
                self.input[0],  # Input FASTA sequence for repeat masking
                seq_db,  # MMseqs database for use in metaeuk or repeat masking
            ]
            if os.path.exists(os.path.join(self.wdir, "tax-report.txt")):
                print("A")
                self.passed_data["tax_assignment"], self.passed_data["tax_id"] = TaxonomyIter.Taxonomy.get_taxonomy(
                    os.path.join(self.wdir, "tax-report.txt"), float(self.cutoff)
                )

        def run(self) -> None:
            super().run()
            # Data to pass forward

        @program_catch
        def run_1(self):
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            seq_db = self.output[1]
            # Create sequence database
            self.log_and_run(
                self.program[
                    "createdb",
                    self.input[0],  # Input FASTA file
                    seq_db,  # Output FASTA sequence db
                ]
            )
            # Run taxonomy search
            self.log_and_run(
                self.program[
                    "taxonomy",
                    seq_db,  # Input FASTA sequence db
                    self.data,  # Input OrthoDB
                    tax_db,  # Output tax db
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags),
                    "--threads", self.threads,
                ]
            )
            # Tax report path
            tax_report = os.path.join(self.wdir, "tax-report.txt")
            # Output results
            self.log_and_run(
                self.program[
                    "taxonomyreport",
                    self.data,  # Input OrthoDB
                    tax_db,  # Input tax db
                    tax_report
                ]
            )
            print("B")
            if self.mode == 0:
                self.passed_data["tax_assignment"], self.passed_data["tax_id"] = "1", "2"
            else:
                self.passed_data["tax_assignment"], self.passed_data["tax_id"] = TaxonomyIter.Taxonomy.get_taxonomy(
                    os.path.join(self.wdir, "tax-report.txt"), float(self.cutoff)
                )

        @staticmethod
        def get_taxonomy(tax_results_file: str, cutoff: float) -> Tuple[str, int]:
            assignment: str = "Eukaryota"  # Default to Eukaryota if nothing better is found
            _id: int = 2759
            if not os.path.exists(tax_results_file):
                return assignment.lower(), _id
            try:
                # Get first line
                _tax_results_file = open(tax_results_file, "r")
                while True:
                    line = next(_tax_results_file).rstrip("\r\n").split("\t")
                    # Parse line for assignment
                    _score, _tax_id, _assignment = float(line[0]), int(line[4]), line[5].lstrip(" ")
                    if _assignment in ("unclassified", "root"):
                        continue
                    if _score >= cutoff:
                        assignment = _assignment
                        _id = _tax_id
                    if _score < cutoff:
                        return assignment.lower(), _id
            except:
                return assignment.lower(), _id
            return assignment.lower(), _id

    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, "taxonomy", *args, **kwargs)


if __name__ == "__main__":
    pass
