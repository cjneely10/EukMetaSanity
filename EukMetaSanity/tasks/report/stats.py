import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch


class ReportStatsIter(TaskList):
    class ReportStats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            for _f in (os.path.join(self.wdir, self.record_id + ".db"),
                       os.path.join(self.wdir, self.record_id + ".summary")):
                if os.path.exists(_f):
                    os.remove(_f)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".db"),
                os.path.join(self.wdir, self.record_id + ".summary"),
            ]

        @program_catch
        def run_1(self):
            # Determine the prefixes to assign to each file based on type
            # Run summary script
            self.log_and_run(
                self.local["summarize_annotations.py"][
                    "-f", self.input[0],
                    "-a", (*self.parse_input()),
                    "-o", os.path.join(self.wdir, self.record_id),
                    "-e", self.max_evalue,
                ]
            )

        def parse_input(self) -> List[str]:
            _out = []
            for _file in self.input:
                if _file.endswith(".emapper"):
                    _out.append("%s=%s" % ("eggnog", _file))
                elif _file.endswith(".kegg"):
                    _out.append("%s=%s" % ("kegg", _file))
                else:
                    _db_name = os.path.splitext(_file)[1][1:-3]
                    if _db_name in ("", "g", "rfam_db"):
                        continue
                    _out.append("%s=%s" % (_db_name, _file))
            return _out

    def __init__(self, *args, **kwargs):
        super().__init__(ReportStatsIter.ReportStats, "stats", *args, **kwargs)


if __name__ == "__main_":
    pass
