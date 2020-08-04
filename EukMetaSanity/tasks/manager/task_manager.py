from typing import Dict, List, Type
from EukMetaSanity.tasks.base.task_class import TaskList
# Run imports
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.run.ab_initio import AbInitioIter
from EukMetaSanity.tasks.run.repeat_modeling import RepeatsIter
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter
# Report imports
from EukMetaSanity.tasks.report.mmseqs import MMseqsIter
from EukMetaSanity.tasks.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.base.summarize import SummarizeIter
from EukMetaSanity.tasks.report.eggnog_mapper import EggNOGMapper
# Refine imports
from EukMetaSanity.tasks.refine.maker import MakerIter
from EukMetaSanity.tasks.refine.assemble import AssembleIter


class TaskManager:
    def __init__(self):
        self._programs = {
            "run": [
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
                EvidenceIter,
                SummarizeIter,
            ],
            "refine": [
                MakerIter,
                AssembleIter,
            ],
            "report": [
                KoFamScanIter,
                MMseqsIter,
                EggNOGMapper,
                SummarizeIter,
            ],
        }
        self._input_type = {
            "refine": ["abinitio", "mask", "nr_gff3", "mask_gff3", "fna", "metaeuk"],
            "report": ["prot"]
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        return self._programs

    @property
    def input_type(self):
        return self._input_type
