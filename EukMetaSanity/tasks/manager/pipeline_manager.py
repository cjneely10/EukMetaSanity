from typing import Dict, List, Type
# General imports
from EukMetaSanity.tasks.base.task_class import TaskList
# Test imports
from EukMetaSanity.tasks.test_tasks.createdb import MMSeqsCreateDBIter
from EukMetaSanity.tasks.test_tasks.searchdb import MMSeqsSearchDBIter

# from EukMetaSanity.tasks.base.summarize import SummarizeIter
# # Run imports
# from EukMetaSanity.tasks.official.run.taxonomy import TaxonomyIter
# from EukMetaSanity.tasks.official.run.ab_initio import AbInitioIter
# from EukMetaSanity.tasks.official.run.repeat_modeling import RepeatsIter
# from EukMetaSanity.tasks.official.run.initial_evidence import EvidenceIter
# # Report imports
# from EukMetaSanity.tasks.official.report.mmseqs import MMseqsIter
# from EukMetaSanity.tasks.official.report.kofamscan import KoFamScanIter
# from EukMetaSanity.tasks.official.report.eggnog_mapper import EggNOGMapper
# from EukMetaSanity.tasks.official.report.stats import ReportStatsIter
# # Refine imports
# from EukMetaSanity.tasks.official.refine.rnaseq import RnaSeqIter
# from EukMetaSanity.tasks.official.refine.braker import BrakerIter
# from EukMetaSanity.tasks.official.refine.transcriptomes import TranscriptomesIter

"""
Class handles ordering of tasks to complete

"""


class PipelineManager:
    def __init__(self):
        self.programs: Dict[str, List[Type[TaskList]]] = {
            # "run": [
            #     TaxonomyIter,
            #     RepeatsIter,
            #     AbInitioIter,
            #     EvidenceIter,
            #     SummarizeIter,
            # ],
            # "report": [
            #     KoFamScanIter,
            #     MMseqsIter,
            #     EggNOGMapper,
            #     ReportStatsIter,
            #     SummarizeIter,
            # ],
            # "refine": [
            #     RnaSeqIter,
            #     TranscriptomesIter,
            #     BrakerIter,
            #     SummarizeIter,
            # ],
            "test": [
                MMSeqsSearchDBIter,
                MMSeqsCreateDBIter,
            ],
        }
