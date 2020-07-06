from typing import Dict, List
from EukMetaSanity.src.tasks.task_class import Task
from EukMetaSanity.src.tasks.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.tasks.tasks.ab_initio import AbInitioIter
from EukMetaSanity.src.tasks.tasks.repeat_modeling import RepeatsIter


class TaskManager:
    def __init__(self):
        self._programs = {
            "run": [
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
            ],
            "refine": [],
            "report": [],
        }

    @property
    def programs(self) -> Dict[str, List[Task]]:
        return self._programs
