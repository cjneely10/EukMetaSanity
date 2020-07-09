import os
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.tasks import Task, TaskList, program_catch

"""
Determine the taxonomy of the Eukaryotic MAG

"""


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            seq_db = os.path.join(self.wdir, self.record_id + "_db")
            results_file = os.path.join(self.wdir, self.record_id + "-tax-report.txt")
            # Expected output
            self.output = {Data.Type.OUT: [
                results_file,  # Taxonomic results for ab initio prediction
                self.input[Data.Type.IN][0],  # Input FASTA file for repeat masking
                seq_db,  # MMseqs database for use in metaeuk or repeat masking
            ]}
            super().run()

        @program_catch
        def run_1(self):
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            seq_db = self.output[Data.Type.OUT][2]
            # Create sequence database
            self.log_and_run(
                self.program[
                    "createdb",
                    self.input[Data.Type.IN][0],  # Input FASTA file
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
            # Output results
            self.log_and_run(
                self.program[
                    "taxonomyreport",
                    self.data,  # Input OrthoDB
                    tax_db,  # Input tax db
                    self.output[Data.Type.OUT][0] + ".tmp"  # Output results file
                ]
            )
            # Write taxonomy result species to file
            with open(self.output[Data.Type.OUT][0], "w") as w:
                w.write(TaxonomyIter.Taxonomy.get_taxonomy(self.output[Data.Type.OUT][0] + ".tmp") + "\n")

        @staticmethod
        def get_taxonomy(tax_results_file: str) -> str:
            assignment: str = "Eukaryota"  # Default to Eukaryota if nothing better is found
            if not os.path.exists(tax_results_file):
                return assignment
            # Get first line
            _tax_results_file = open(tax_results_file, "r")
            try:
                while True:
                    line = next(_tax_results_file).rstrip("\r\n").split("\t")
                    # Parse line for assignment
                    _score, _tax_id, _assignment = float(line[0]), line[4], line[5].replace(" ", "")
                    if _assignment in ("unclassified", "root"):
                        continue
                    if _score < 80.0:
                        break
                    # Keep new value if >= 80.0% of contigs map to the taxonomy
                    else:
                        assignment = _assignment
            except StopIteration:
                return assignment

    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, "taxonomy", *args, **kwargs)


if __name__ == "__main__":
    pass
