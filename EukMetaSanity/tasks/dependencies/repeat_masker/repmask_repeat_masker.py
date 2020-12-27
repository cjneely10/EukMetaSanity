import os
from pathlib import Path
from EukMetaSanity import ProcessExecutionError, set_complete
from EukMetaSanity import Task, TaskList, program_catch, prefix, DependencyInput


class RepeatMaskerIter(TaskList):
    name = "repmask.repeat_masker"
    requires = ["taxonomy"]
    depends = [DependencyInput("repmod.repeat_modeler")]
    
    class RepeatMasker(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            data_files = []
            data_files += [_f for _f in self.data if _f != ""]
            # Perform on optimal taxonomic identification
            if self.input["taxonomy"]["taxonomy"].family is not None:
                data_files += [self.input["taxonomy"]["taxonomy"].family.value]
            _file = str(self.input["repmod.repeat_modeler"]["model"])
            if os.path.exists(_file) and os.path.getsize(_file) > 0:
                data_files.append(_file)
            out = []
            for _search in data_files:
                if "classified" in _search:
                    search = ("-lib", str(Path(_search).resolve()))
                    _dir = os.path.join(self.wdir, "repeats_" + prefix(_search))
                else:
                    search = ("-species", _search)
                    _dir = os.path.join(self.wdir, "repeats_" + _search.replace(" ", "_"))
                out.append((search, _dir))
            self.output = {
                "libraries": out
            }
            
        @program_catch
        def run(self):
            _added_dirs = []
            for val in self.output["libraries"]:
                search, _dir = val
                # Create contained directory
                if os.path.exists(_dir):
                    continue
                os.makedirs(_dir)
                # Call RepeatMasker on modeled repeats in the new directory
                script = self.create_script(
                    self.program[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", _dir,
                        self.input["root"]["fna"],
                    ],
                    "%s.sh" % os.path.basename(_dir)
                )
                try:
                    self.parallel(script)
                except ProcessExecutionError:
                    continue
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerIter.RepeatMasker, RepeatMaskerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
