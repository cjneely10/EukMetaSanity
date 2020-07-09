from typing import Dict, List
from EukMetaSanity.src.tasks.base.task_class import Task
from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.tasks.ab_initio import AbInitioIter
from EukMetaSanity.src.tasks.repeat_modeling import RepeatsIter
from EukMetaSanity.src.tasks.initial_evidence import EvidenceIter


class TaskManager:
    def __init__(self):
        self._programs = {
            "run": [
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
                EvidenceIter,
            ],
            "refine": [],
            "report": [],
        }

    @property
    def programs(self) -> Dict[str, List[Task]]:
        return self._programs
