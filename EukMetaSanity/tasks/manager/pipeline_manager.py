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

"""
Class handles ordering of tasks to complete

"""


class PipelineManager:
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
            ]
        }
