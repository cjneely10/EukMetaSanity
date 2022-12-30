import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from EukMetaSanity.testing_imports import MMSeqsTaxonomyReportParser

single_branch = """13.5938\t755\t755\tno rank\t0\tunclassified
86.4062\t4799\t0\tno rank\t1\troot
86.3702\t4797\t0\tno rank\t131567\t  cellular organisms
86.1001\t4782\t25\tsuperkingdom\t2759\t    Eukaryota
"""

double_branch = single_branch + """5.5\t100\t4\tno rank\t131568\t  uncellular organisms
1.1\t25\t1\tsuperkingdom\t2759\t    Uncellulota
"""


class TestTaxonomyReportParser(unittest.TestCase):
    def test_single_branch(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(single_branch)
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
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(double_branch)
        tree = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(tree.children[1].children) == 2
        assert tree.children[1].children[1].data.scientific_name == "uncellular organisms"
        assert len(tree.children[1].children[1].children) == 1
        assert tree.children[1].children[1].children[0].data.scientific_name == "Uncellulota"

    def test_find_best_taxonomy(self):
        pass
