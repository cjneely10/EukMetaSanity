"""
Module describe outer-most Tasks to complete for a given pipeline
"""

from typing import Dict, List, Type
# General imports
from EukMetaSanity.tasks.base.task_class import TaskList
# # Run imports
from EukMetaSanity.tasks.official.report.quality import QualityIter
from EukMetaSanity.tasks.official.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.official.run.repeats import RepeatsIter
from EukMetaSanity.tasks.official.run.abinitio_augustus import AbInitioAugustusIter
from EukMetaSanity.tasks.official.run.abinitio_genemark import AbInitioGeneMarkIter
from EukMetaSanity.tasks.official.run.merge import MergeIter
from EukMetaSanity.tasks.official.run.evidence import EvidenceIter
# # Report imports
from EukMetaSanity.tasks.official.report.kofamscan import KoFamScanIter
from EukMetaSanity.tasks.official.report.eggnog import EggNogMapperIter
from EukMetaSanity.tasks.official.report.mmseqs import MMseqsIter
# # Refine imports
from EukMetaSanity.tasks.official.refine.mapping import MappingIter
from EukMetaSanity.tasks.official.refine.filtering import FilteringIter


# pylint: disable=too-few-public-methods
class PipelineManager:
    """ Class handles ordering of tasks to complete

    Key represents the program name that the user will pass on the command line

    Example:
    self.program["test"] = [MMSeqsSearchDBIter, MMSeqsCreateDBIter]

    """
    def __init__(self):
        """ Load available pipeline data

        """
        self._programs: Dict[str, List[Type[TaskList]]] = {
            "run": [
                MergeIter,
                TaxonomyIter,
                AbInitioAugustusIter,
                AbInitioGeneMarkIter,
                RepeatsIter,
                EvidenceIter,
            ],
            "report": [
                KoFamScanIter,
                EggNogMapperIter,
                MMseqsIter,
                QualityIter,
            ],
            "refine": [
                MappingIter,
                FilteringIter
            ]
            # Add any user-created pipelines here!
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        """
        Get pipeline programs available
        """
        return self._programs
