from typing import Dict, List, Type
# General imports
from EukMetaSanity.tasks.base.task_class import TaskList
# Test imports
from EukMetaSanity.tasks.test_tasks.createdb import MMSeqsCreateDBIter
from EukMetaSanity.tasks.test_tasks.searchdb import MMSeqsSearchDBIter
# # # Test new api
# # Run imports
from EukMetaSanity.tasks.test_tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.test_tasks.run.repeats import RepeatsIter
from EukMetaSanity.tasks.test_tasks.run.abinitio import AbInitioIter
from EukMetaSanity.tasks.test_tasks.run.evidence import EvidenceIter
# # Report imports
from EukMetaSanity.tasks.test_tasks.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.test_tasks.report.eggnog import EggNogMapperIter
from EukMetaSanity.tasks.test_tasks.report.mmseqs import MMseqsIter
from EukMetaSanity.tasks.test_tasks.report.stats import StatsIter


class PipelineManager:
    """ Class handles ordering of tasks to complete

    Key represents the program name that the user will pass on the command line

    :param programs: Programs in pipeline (need not be in any order)
    :type programs: Dict[str, List[Type[TaskList]]]
    """
    def __init__(self):
        self.programs: Dict[str, List[Type[TaskList]]] = {
            "test": [
                MMSeqsSearchDBIter,
                MMSeqsCreateDBIter,
            ],
            "test_run": [
                EvidenceIter,
                TaxonomyIter,
                AbInitioIter,
                RepeatsIter,
            ],
            "test_report": [
                StatsIter,
                KoFamScanIter,
                EggNogMapperIter,
                MMseqsIter,
            ],
        }
