from typing import Dict, List
from EukMetaSanity.tasks.base.task_class import Task
from EukMetaSanity.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.ab_initio import AbInitioIter
from EukMetaSanity.tasks.add_results import AddResultsIter
from EukMetaSanity.tasks.repeat_modeling import RepeatsIter
from EukMetaSanity.tasks.initial_evidence import EvidenceIter


class TaskManager:
    def __init__(self):
        self._programs = {
            "run": [
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
                EvidenceIter,
                AddResultsIter,
            ],
            "refine": [],
            "report": [],
        }

    @property
    def programs(self) -> Dict[str, List[Task]]:
        return self._programs
