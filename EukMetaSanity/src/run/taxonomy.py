import os
from pathlib import Path
from typing import List, Union, Type

from yapim import Task, DependencyInput

from .mmseqs_taxonomy_report_parser import MMSeqsTaxonomyReportParser


class Taxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": self.input["MMSeqsTaxonomy"]["tax-report"],
            "final": ["tax-report", "taxonomy"]
        }
        self.set_tax_file()

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return ["MetaEukEV"]

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsTaxonomy", {"MetaEukEV": {"evidence-prot": "fasta"}})]

    def run(self):
        self.set_tax_file()

    def set_tax_file(self):
        tax_file = self.output["tax-report"]
        if os.path.exists(tax_file):
            self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(Path(tax_file).resolve())
        else:
            self.output["taxonomy"] = []
