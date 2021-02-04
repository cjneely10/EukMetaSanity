"""
Module holds logic to identify taxonomy of a genome/MAG
"""
from EukMetaSanity import TaxonomyAssignment
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, set_complete


class TaxonomyIter(TaskList):
    """ This class will use `mmseqs` to identify putative taxonomy for an organism.

    name: taxonomy

    requires:

    depends: mmseqs.taxonomy

    output: taxonomy[TaxonomyAssignment]

    final: mmseqs.taxonomy.tax-report[Path], taxonomy[TaxonomyAssignment]

    config:
        cutoff: 8.0  # Minimum % of mapped reads to tax level

    """
    name = "taxonomy"
    requires = []
    depends = [DependencyInput("mmseqs.taxonomy")]

    class Taxonomy(Task):
        """
        Predict taxonomy
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "taxonomy": self.get_taxonomy(
                    str(self.input["mmseqs.taxonomy"]["tax-report"]),
                    float(self.config["cutoff"])),
                "final": ["mmseqs.taxonomy.tax-report", "taxonomy"]
            }

        @program_catch
        def run(self):
            """
            Run
            """
            pass

        def get_taxonomy(self, tax_results_file: str, cutoff: float) -> TaxonomyAssignment:
            """ Parse taxonomy results into TaxonomyAssignment object for downstream use

            :param tax_results_file: File containing MMseqs taxonomy output summary results
            :param cutoff: Minimum percent of reads to allow for parsing a taxonomic level
            :return: TaxonomyAssignment object containing results
            """
            if self.developer_mode:
                return TaxonomyAssignment()
            tax_assignment_out = TaxonomyAssignment()
            tax_levels = ["kingdom", "phylum", "class", "order", "superfamily", "family", "genus", "species"]
            taxonomy = {key: None for key in tax_levels}
            tax_assignment_out.kingdom = TaxonomyAssignment.Assignment("Eukaryota", 2759, -1)
            try:
                # Get first line
                _tax_results_file = open(tax_results_file, "r")
                _level = ""
                while True:
                    line = next(_tax_results_file).rstrip("\r\n").split("\t")
                    # Parse line for assignment
                    _score, _level, _tax_id, _assignment = float(line[0]), line[3], int(line[4]), line[5].lstrip(" ")
                    if _assignment in ("unclassified", "root"):
                        continue
                    if _score >= cutoff and taxonomy.get(_level, None) is None:
                        taxonomy[_level] = TaxonomyAssignment.Assignment(_assignment, _tax_id, _score)
            except StopIteration:
                pass
            for level, assignment in taxonomy.items():
                if level == "class":
                    setattr(tax_assignment_out, "_class", assignment)
                else:
                    setattr(tax_assignment_out, level, assignment)
            return tax_assignment_out

    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)
