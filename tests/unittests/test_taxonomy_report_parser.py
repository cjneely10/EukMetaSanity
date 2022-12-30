import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from EukMetaSanity.testing_imports import MMSeqsTaxonomyReportParser


class TestTaxonomyReportParser(unittest.TestCase):
    file = Path(__file__).parent.joinpath("data").joinpath("test-taxonomy.txt")

    def test_single_branch(self):
        temp_file = NamedTemporaryFile()
        with open(temp_file.name, "w") as temp_file_ptr:
            temp_file_ptr.write("""13.5938	755	755	no rank	0	unclassified
86.4062	4799	0	no rank	1	root
86.3702	4797	0	no rank	131567	  cellular organisms
86.1001	4782	25	superkingdom	2759	    Eukaryota
77.6737	4314	3	clade	33154	      Opisthokonta
77.0436	4279	1	kingdom	33208	        Metazoa
76.9175	4272	7	clade	6072	          Eumetazoa
75.9093	4216	25	clade	33213	            Bilateria
66.8347	3712	3	clade	33317	              Protostomia
65.1242	3617	2	clade	1206794	                Ecdysozoa
64.6201	3589	0	clade	88770	                  Panarthropoda
64.6201	3589	11	phylum	6656	                    Arthropoda
60.8570	3380	2	clade	197563	                      Mandibulata
60.5149	3361	8	clade	197562	                        Pancrustacea
""")
        tree = MMSeqsTaxonomyReportParser.create_tree(Path(temp_file.name))
        assert tree.children[0].data.scientific_name_cleaned == "unclassified"
        assert tree.children[1].data.scientific_name_cleaned == "root"
        assert tree.children[1].children[0].data.scientific_name_cleaned == "cellular organisms"
        assert tree.children[1].children[0].children[0].data.scientific_name_cleaned == "Eukaryota"

    def test_empty(self):
        pass

    def test_two_branches(self):
        pass
