import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch


class StatsIter(TaskList):
    name = "stats"
    requires = ["emapper", "kofamscan", "mmseqs"]
    
    class Stats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            for _f in (os.path.join(self.wdir, self.record_id + ".db"),
                       os.path.join(self.wdir, self.record_id + ".summary")):
                if os.path.exists(_f):
                    os.remove(_f)
            self.output = {
                "summary-db": os.path.join(self.wdir, self.record_id + ".db"),
                "summary": os.path.join(self.wdir, self.record_id + ".summary"),
            }
            
        @program_catch
        def run(self):
            # Determine the prefixes to assign to each file based on type
            # Run summary script
            self.local["summarize_annotations.py"][
                "-f", str(self.input["root"]["prot"]),
                "-a", (*self.parse_input()),
                "-o", os.path.join(self.wdir, self.record_id),
                "-e", self.max_evalue,
            ]()

        def _parse_input(self) -> List[str]:
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__(StatsIter.Stats, StatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
