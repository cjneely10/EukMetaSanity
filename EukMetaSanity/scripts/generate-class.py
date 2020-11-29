#!/usr/bin/env python3
import os
from pathlib import Path
from EukMetaSanity.utils.arg_parse import ArgParse


boilerplate = """import os
from EukMetaSanity import Task, TaskList, program_catch, prefix, touch


class {0}Iter(TaskList):
    name = "{1}"
    requires = []
    
    class {0}(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {2}
            
        @program_catch
        def run(self):
            pass
            
    def __init__(self, *args, **kwargs):
        super().__init__({0}Iter.{0}, {0}Iter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
"""

config_boilerplate = """[SLURM]
# Set to True if using SLURM
USE_CLUSTER = False
## WORKERS on each config section will assign the number of jobs to run/MAGs to analyze simultaneously
## THREADS on each config section set the number of cores to run for each job
## MEMORY on each config section sets memory size requests. Make sure this matches mmseqs requested memory
## Pass any flags you wish below
## DO NOT PASS the following: --nodes, --ntasks, --mem
--qos = unlim
--job-name = EukMS
user-id = uid

"""


def create_class_file(file_path: str, class_name: str, cfg_name: str):
    fp = open(os.path.join(file_path, cfg_name + ".py"), "w")
    fp.write(boilerplate.format(class_name, cfg_name, """{
                
            }"""))
    fp.close()


def update_config_file(cfg_path: str, cfg_name: str):
    cfg_path = str(Path(cfg_path).resolve())
    if not os.path.exists(cfg_path):
        fp = open(cfg_path, "w")
        fp.write(config_boilerplate)
        fp.close()
    fp = open(cfg_path, "a")
    fp.write("[%s]\nWORKERS = 1\nTHREADS = 1\nMEMORY = 10G\nTIME = 1:00:00\nFLAGS = ''\n\n" % cfg_name)
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
            (("-c", "--config_path"),
             {"help": "Update existing base config file to include new class, or create fresh config file"}),
        ),
        description="Create new EukMetaSanity class"
    )
    create_class_file(ap.args.path, ap.args.class_name, ap.args.config_name)
    if ap.args.config_path is not None:
        update_config_file(ap.args.config_path, ap.args.config_name)
