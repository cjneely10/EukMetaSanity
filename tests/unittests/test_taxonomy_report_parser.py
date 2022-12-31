import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from EukMetaSanity.testing_imports import MMSeqsTaxonomyReportParser

unclassified = "100.0\t1000\t1000\tno rank\t0\tunclassified\n"

single_branch = """13.5938\t755\t755\tno rank\t0\tunclassified
86.4062\t4799\t0\tno rank\t1\troot
86.3702\t4797\t0\tno rank\t131567\t  cellular organisms
86.1001\t4782\t25\tsuperkingdom\t2759\t    Eukaryota
"""

nested = """18.1329\t202\t202\tno rank\t0\tunclassified
81.8671\t912\t0\tno rank\t1\troot
81.7774\t911\t0\tno rank\t131567\t  cellular organisms
81.5081\t908\t6\tsuperkingdom\t2759\t    Eukaryota
73.8779\t823\t1\tclade\t33154\t      Opisthokonta
72.5314\t808\t0\tkingdom\t33208\t        Metazoa
72.3519\t806\t3\tclade\t6072\t          Eumetazoa
81.5081\t908\t4\tsuperkingdom\t2759\t    Archaea
73.8779\t823\t2\tclade\t33154\t      Archaea-clade
72.5314\t808\t0\tkingdom\t33208\t        Archaea-kingdom
72.3519\t806\t5\tclade\t6072\t          Archaea-clade2
"""


class TestTaxonomyReportParser(unittest.TestCase):
    def test_single_branch(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(single_branch)
        tree, leaf_nodes = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(leaf_nodes) == 2
        assert len(tree.children) == 2
        assert tree.children[0].data.scientific_name == "unclassified"
        assert tree.children[1].data.scientific_name == "root"
        assert len(tree.children[1].children) == 1
        assert tree.children[1].children[0].data.scientific_name == "cellular organisms"
        assert len(tree.children[1].children[0].children) == 1
        assert tree.children[1].children[0].children[0].data.scientific_name == "Eukaryota"

        assert tree.children[0].cumulative_read_count == 755
        assert tree.children[1].cumulative_read_count == 0
        assert tree.children[1].children[0].cumulative_read_count == 0
        assert tree.children[1].children[0].children[0].cumulative_read_count == 25

        assert MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root", "cellular organisms", "Eukaryota")
        assert not MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root", "cellular organisms", "Eukaryota", "Bil")

    def test_empty(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write("")
        tree, leaf_nodes = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(leaf_nodes) == 0
        assert len(tree.children) == 0
        assert not MMSeqsTaxonomyReportParser._find_taxonomy(tree, "root")

    def test_unclassified(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(unclassified)
        tree, leaf_nodes = MMSeqsTaxonomyReportParser._create_tree(Path(temp_file.name))
        assert len(leaf_nodes) == 1
        assert len(tree.children) == 1
        assert tree.children[0].data.scientific_name == "unclassified"

    def test_find_best_taxonomy(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write(nested)
        best_taxonomy = MMSeqsTaxonomyReportParser.find_best_taxonomy(Path(temp_file.name))
        assert best_taxonomy[-1][1]["value"] == "Archaea-clade2"

        assigned_kingdom = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(best_taxonomy, "kingdom")
        assert assigned_kingdom[0] == "kingdom"
        assert assigned_kingdom[1]["value"] == "Archaea-kingdom"

        request_for_class = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(best_taxonomy, "class")
        assert request_for_class[0] == "kingdom"
        assert request_for_class[1]["value"] == "Archaea-kingdom"

    def test_real_file(self):
        file = Path(__file__).parent.joinpath("data").joinpath("test-taxonomy.txt")
        best_taxonomy = MMSeqsTaxonomyReportParser.find_best_taxonomy(file)
        assert best_taxonomy[-1][1]["value"] == "Tigriopus californicus"

        assigned_phylum = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(best_taxonomy, "phylum")
        assert assigned_phylum[0] == "phylum"
        assert assigned_phylum[1]["value"] == "Arthropoda"

        request_for_sf = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(best_taxonomy, "superfamily")
        assert request_for_sf[0] == "order"
        assert request_for_sf[1]["value"] == "Harpacticoida"

        invalid_search = MMSeqsTaxonomyReportParser.find_assignment_nearest_request(best_taxonomy, "meow")
        assert invalid_search[0] == "superkingdom"
        assert invalid_search[1]["value"] == "Eukaryota"
