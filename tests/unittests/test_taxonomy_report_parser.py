import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from EukMetaSanity.testing_imports import MMSeqsTaxonomyReportParser


class TestTaxonomyReportParser(unittest.TestCase):
    file = Path(__file__).parent.joinpath("data").joinpath("test-taxonomy.txt")
    single_branch = """13.5938	755	755	no rank	0	unclassified
86.4062	4799	0	no rank	1	root
86.3702	4797	0	no rank	131567	  cellular organisms
86.1001	4782	25	superkingdom	2759	    Eukaryota
"""

    def test_single_branch(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(TestTaxonomyReportParser.single_branch)
        tree = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(tree.children) == 2
        assert tree.children[0].data.scientific_name == "unclassified"
        assert tree.children[1].data.scientific_name == "root"
        assert len(tree.children[1].children) == 1
        assert tree.children[1].children[0].data.scientific_name == "cellular organisms"
        assert len(tree.children[1].children[0].children) == 1
        assert tree.children[1].children[0].children[0].data.scientific_name == "Eukaryota"

        assert MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root", "cellular organisms", "Eukaryota")
        assert not MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root", "cellular organisms", "Eukaryota", "Bil")

    def test_empty(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write("")
        tree = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(tree.children) == 0
        assert not MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root")

    def test_two_branches(self):
        pass

    def test_find_in_file(self):
        pass

    def test_find_best_taxonomy(self):
        pass
