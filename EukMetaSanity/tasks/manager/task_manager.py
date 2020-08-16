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
# Fast_refine imports
from EukMetaSanity.tasks.fast_refine.rnaseq import RnaSeqIter
from EukMetaSanity.tasks.fast_refine.transcriptomes import TranscriptomesIter
from EukMetaSanity.tasks.fast_refine.merge import MergeIter

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
            "fast_refine": [
                RnaSeqIter,
                TranscriptomesIter,
                MergeIter,
                SummarizeIter,
            ],
        }
        self._input_type = {
            "report": ["prot", "nr_gff3", "fna"],
            "fast_refine": ["fna", "nr_gff3"],
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        return self._programs

    @property
    def input_type(self):
        return self._input_type
