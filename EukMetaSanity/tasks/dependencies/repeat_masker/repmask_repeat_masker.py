import os
from pathlib import Path
from EukMetaSanity import ProcessExecutionError
from EukMetaSanity import Task, TaskList, program_catch, prefix


class RepeatMaskerIter(TaskList):
    name = "repmask.repeat_masker"
    requires = ["taxonomy"]
    depends = ["repmod.repeat_modeler"]
    
    class RepeatMasker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            _file = str(self.input["repmod.repeat_modeler"]["model"])
            data_files = []
            data_files += [_f for _f in self.data if _f != ""]
            # Perform on optimal taxonomic identification
            if self.input["taxonomy"]["taxonomy"].family is not None:
                data_files += [self.input["taxonomy"]["taxonomy"].family.value]
            if os.path.exists(_file) and os.path.getsize(_file) > 0:
                data_files.append(_file)
            self.output = {
                "libraries": data_files
            }
            
        @program_catch
        def run(self):
            _added_dirs = []
            for _search in self.output["libraries"]:
                # Parse for if as file or a RepeatMasker library
                if "RM" in _search:
                    search = ("-lib", str(Path(_search).resolve()))
                    _dir = "repeats_" + prefix(_search)
                else:
                    search = ("-species", _search)
                    _dir = "repeats_" + _search.replace(" ", "_")
                # Create contained directory
                self.pm.add_dirs(self.record_id, [_dir])
                _added_dirs.append(self.pm.get_dir(self.record_id, _dir))
                # Call RepeatMasker on modeled repeats in the new directory
                script = self.create_script(
                    self.program[
                        "-pa", self.threads,
                        (*self.added_flags),
                        (*search),
                        "-dir", self.pm.get_dir(self.record_id, _dir),
                        self.input["root"]["fna"],
                    ],
                    "mask.sh"
                )
                try:
                    self.parallel(script)
                except ProcessExecutionError:
                    continue
            
    def __init__(self, *args, **kwargs):
        super().__init__(RepeatMaskerIter.RepeatMasker, RepeatMaskerIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
