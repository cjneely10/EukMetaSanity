from typing import Dict, List, Type
# General imports
from EukMetaSanity.tasks.base.task_class import TaskList
# # Run imports
from EukMetaSanity.tasks.official.run.taxonomy import TaxonomyIter
# from EukMetaSanity.tasks.official.run.repeats import RepeatsIter
# from EukMetaSanity.tasks.official.run.abinitio import AbInitioIter
# from EukMetaSanity.tasks.official.run.evidence import EvidenceIter
# # Report imports
from EukMetaSanity.tasks.test_tasks.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.test_tasks.report.eggnog import EggNogMapperIter
from EukMetaSanity.tasks.test_tasks.report.mmseqs import MMseqsIter
from EukMetaSanity.tasks.test_tasks.report.stats import StatsIter


class PipelineManager:
    """ Class handles ordering of tasks to complete

    Key represents the program name that the user will pass on the command line

    Example:
    self.program["test"] = [MMSeqsSearchDBIter, MMSeqsCreateDBIter]

    """
    def __init__(self):
        self.programs: Dict[str, List[Type[TaskList]]] = {
            "run": [
                # EvidenceIter,
                TaxonomyIter,
                # AbInitioIter,
                # RepeatsIter,
            ],
            "test_report": [
                StatsIter,
                KoFamScanIter,
                EggNogMapperIter,
                MMseqsIter,
            ],
        }
