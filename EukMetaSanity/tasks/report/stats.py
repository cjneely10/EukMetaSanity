import os
from typing import List
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class ReportStatsIter(TaskList):
    class ReportStats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".db"),
                os.path.join(self.wdir, self.record_id + ".summary"),
            ]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            # Determine the prefixes to assign to each file based on type
            # Run summary script
            self.log_and_run(
                self.local["create-final-annotations.py"][
                    "-f", self.input[0],
                    "-a", *self.parse_input(),
                ]
            )

        def parse_input(self) -> List[str]:
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(ReportStatsIter.ReportStats, "stats", *args, **kwargs)


if __name__ == "__main_":
    pass
