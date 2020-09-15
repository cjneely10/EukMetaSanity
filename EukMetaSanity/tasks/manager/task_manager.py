from typing import Dict, List, Type
from EukMetaSanity.tasks.base.task_class import TaskList
# General imports
from EukMetaSanity.tasks.base.summarize import SummarizeIter
# Run imports
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.run.ab_initio import AbInitioIter
from EukMetaSanity.tasks.run.repeat_modeling import RepeatsIter
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter
# Report imports
from EukMetaSanity.tasks.report.mmseqs import MMseqsIter
from EukMetaSanity.tasks.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.report.eggnog_mapper import EggNOGMapper
from EukMetaSanity.tasks.report.stats import ReportStatsIter
# Refine imports
from EukMetaSanity.tasks.refine.rnaseq import RnaSeqIter
from EukMetaSanity.tasks.refine.braker import BrakerIter
from EukMetaSanity.tasks.refine.transcriptomes import TranscriptomesIter

"""
Class handles ordering of tasks to complete

"""


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
            "report": [
                KoFamScanIter,
                MMseqsIter,
                EggNOGMapper,
                ReportStatsIter,
                SummarizeIter,
            ],
            "refine": [
                RnaSeqIter,
                TranscriptomesIter,
                BrakerIter,
                SummarizeIter,
            ],
        }
        self._input_type = {
            "report": ["prot", "nr_gff3", "fna"],
            "refine": ["mask", "tax", "fna", "abinitio", "metaeuk"],
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        return self._programs

    @property
    def input_type(self):
        return self._input_type
