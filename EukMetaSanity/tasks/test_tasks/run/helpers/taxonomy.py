from typing import Tuple


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
