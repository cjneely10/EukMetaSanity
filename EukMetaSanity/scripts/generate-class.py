#!/usr/bin/env python3
import os
from pathlib import Path
from EukMetaSanity.utils.arg_parse import ArgParse


boilerplate = """import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix, touch


class {0}Iter(TaskList):
    class {0}(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            
        @program_catch
        def run_1(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__({0}Iter.{0}, "{1}", [], 
                         *args, **kwargs)


if __name__ == "__main_":
    pass
"""


def create_file(file_path: str, name: str, cfg_name: str):
    fp = open(os.path.join(file_path, cfg_name + ".py"), "w")
    fp.write(boilerplate.format(name, cfg_name))
    fp.close()


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("class_name",),
             {"help": "Name of class to create"}),
            (("config_name",),
             {"help": "Matching config file section name"}),
            (("-p", "--path"),
             {"help": "Location to generate, default current directory", "default": os.getcwd()}),
        ),
        description="Create new EukMetaSanity class"
    )
    create_file(ap.args.path, ap.args.class_name, ap.args.config_name)
