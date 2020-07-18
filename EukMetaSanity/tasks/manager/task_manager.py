from typing import Dict, List, Type
from EukMetaSanity.tasks.base.task_class import TaskList
# Base summarize function
from EukMetaSanity.tasks.manager.summarize import SummarizeIter
# Run imports
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.run.ab_initio import AbInitioIter
from EukMetaSanity.tasks.run.repeat_modeling import RepeatsIter
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter
# Report imports
from EukMetaSanity.tasks.report.pfam import PfamIter
from EukMetaSanity.tasks.report.mmseqs import MMseqsIter
from EukMetaSanity.tasks.report.hhsuite import HHsuiteIter
from EukMetaSanity.tasks.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.report.eggnog_mapper import EggNOGMapper


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
            "refine": [],
            "report": [
                KoFamScanIter,
                MMseqsIter,
                HHsuiteIter,
                PfamIter,
                EggNOGMapper,
                SummarizeIter,
            ],
        }
        self._input_type = {
            "run": None,
            "refine": "all",  # TODO: Handle code for all step
            "report": "prot"
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        return self._programs

    @property
    def input_type(self):
        return self._input_type
