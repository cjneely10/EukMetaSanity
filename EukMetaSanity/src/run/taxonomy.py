from pathlib import Path
from typing import List, Union, Type

from yapim import Task, DependencyInput, touch

from EukMetaSanity.mmseqs_taxonomy_report_parser import MMSeqsTaxonomyReportParser


class Taxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": self.input["MMSeqsTaxonomy"]["tax-report"],
            "final": ["tax-report", "taxonomy"],
            "_": self.wdir.joinpath(".done")
        }
        tax_file = Path(self.output["tax-report"]).resolve()
        if tax_file.exists():
            self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(tax_file)

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["MetaEukEV", "CreateMMSeqsDB"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsTaxonomy",
                                {"MetaEukEV": {"evidence-prot": "fasta"}, "CreateMMSeqsDB": ["db"]})]

    def run(self):
        tax_file = Path(self.output["tax-report"]).resolve()
        self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(tax_file)
        touch(str(self.output["_"]))
