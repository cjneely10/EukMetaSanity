from collections import namedtuple
from typing import Tuple, Optional
from EukMetaSanity import Task, TaskList, program_catch

Assignment = namedtuple("Assignment", ("value", "score"))


class TaxonomyAssignment:
    """ This class holds the results of a taxonomic assignment by varying taxonomic levels

    """
    def __init__(self):
        self.kingdom: Optional[Assignment] = None
        self.phylum: Optional[Assignment] = None
        self._class: Optional[Assignment] = None
        self.order: Optional[Assignment] = None
        self.superfamily: Optional[Assignment] = None
        self.family: Optional[Assignment] = None
        self.genus: Optional[Assignment] = None
        self.species: Optional[Assignment] = None

    def assignment(self, level: str) -> Optional[Assignment]:
        """ Get Assignment object at given level.
        Returns None is level not found in file, or if level did not exist in file

        :param level: Taxonomy level string (e.g. kingdom, order, _class, etc.)
        :return: Assignment object or None
        """
        return getattr(self, level, None)


class TaxonomyIter(TaskList):
    """ This class will use `mmseqs` to identify putative taxonomy for an organism.

    name: taxonomy

    requires: mmseqs.createdb, mmseqs.taxonomy

    output keys: taxonomy

    finalizes: mmseqs.taxonomy.tax-report

    """
    name = "taxonomy"
    requires = ["mmseqs.createdb", "mmseqs.taxonomy"]
    
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "taxonomy": TaxonomyAssignment(),
                "final": ["mmseqs.taxonomy.tax-report"]
            }
            
        @program_catch
        def run(self):
            pass

        @staticmethod
        def get_taxonomy(tax_results_file: str, cutoff: float, deepest_level: str = "strain") -> TaxonomyAssignment:
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
