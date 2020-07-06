from EukMetaSanity.src.tasks.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.tasks.tasks.repeat_modeling import RepeatsIter
from EukMetaSanity.src.tasks.tasks.ab_initio import AbInitioIter


class TaskManager:
    def __init__(self):
        self.tasks = {
            "run": (
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
            ),
            "refine": (),
            "report": (),
        }
