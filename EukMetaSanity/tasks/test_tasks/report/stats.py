import os
from typing import List
from EukMetaSanity import Task, TaskList, program_catch


class StatsIter(TaskList):
    """ Task summarizes all gene call annotations into a tsv and sqlite3 database

    Outputs: summary-db, summary
    Finalizes: summary-db, summary

    """
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
                "final": ["summary-db", "summary"]
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

        def parse_input(self) -> List[str]:
            _out = []
            for _file in (
                str(self.input["emapper"]["emapper"]),
                str(self.input["kofamscan"]["kegg"]),
                *(str(v) for k, v in self.input["mmseqs"].items() if k != "final"),
            ):
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
        super().__init__(StatsIter.Stats, StatsIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
