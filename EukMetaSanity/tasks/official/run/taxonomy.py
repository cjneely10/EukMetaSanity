from EukMetaSanity import Task, TaskList, program_catch


class TaxonomyIter(TaskList):
    """ This class will use `mmseqs` to identify putative taxonomy for an organism.

    Outputs: seq_db, tax_db, tax-report

    Finalizes: tax-report

    """
    name = "taxonomy"
    requires = ["mmseqs.createdb", "mmseqs.taxonomy"]
    
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "Parser": TaxonomyIter.Taxonomy.get_taxonomy,
                "seq_db": self.input["mmseqs.createdb"]["db"],
                "tax-report": self.input["mmseqs.taxonomy"]["tax-report"],
                "final": ["tax-report"]
            }
            
        @program_catch
        def run(self):
            pass

        @staticmethod
        def get_taxonomy(tax_results_file: str, cutoff: float, deepest_level: str = "strain") -> Tuple[str, int]:
            tax_levels = ["kingdom", "phylum", "class", "order", "superfamily", "family", "genus", "species"]
            taxonomy = {key: None for key in tax_levels}
            assert deepest_level in taxonomy.keys()
            assignment: str = "Eukaryota"  # Default to Eukaryota if nothing better is found
            _id: int = 2759
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
                        taxonomy[_level] = (_assignment, _tax_id)
            except StopIteration:
                pass
            # Found requested level
            if taxonomy.get(deepest_level, None) is not None:
                return taxonomy[deepest_level]
            # Return deepest level with matching cutoff
            for i in range(len(tax_levels) - 1, -1, -1):
                if taxonomy[tax_levels[i]] is not None:
                    return taxonomy[tax_levels[i]]
            return '"%s"' % assignment.lower(), _id
            
    def __init__(self, *args, **kwargs):
        super().__init__(TaxonomyIter.Taxonomy, TaxonomyIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
