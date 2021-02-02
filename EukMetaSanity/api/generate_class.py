#!/usr/bin/env python

"""
Module holds logic to generate TaskList/Task class stub and associated config file sections
"""

import os
from typing import Optional, List
import yaml

from EukMetaSanity.arg_parse import ArgParse


BOILERPLATE = '''"""
Module holds {1} build functionality
"""

import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class {0}Iter(TaskList):
    """ TaskList class iterates over {1} tasks

    name: {1}

    requires:

    depends:

    expects:

    output:

    final:

    config:

    """
    name = "{1}"
    requires = []
    depends = []

    class {0}(Task):
        """
        Task class handles {1} task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {2}

        @program_catch
        def run(self):
            """
            Run {1}
            """
            pass

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__({0}Iter.{0}, {0}Iter.name, *args, **kwargs)
'''

CONFIG_OUTER_LEVEL = {
    "workers": 1,
    "threads": 1,
    "memory": "16G",
    "time": "1:00:00",
    "dependencies": {}
}

CONFIG_INNER_LEVEL = {
    "program": "",
    "data": "",
    "FLAGS": []
}

NEW_CONFIG_DATA = {
    "INPUT": {
        "base": "root"
    },
    "SLURM": {
        "USE_CLUSTER": "false",
        "--qos": "unlim",
        "--job-name": "EukMS",
        "user-id": "uid"
    }
}


def update_config_file(file_path: str, cfg_name: str, existing_sections: Optional[List[str]] = None):
    """ Add the newly-created class base config section at either root or within existing section

    :param file_path: Config file path
    :param cfg_name: Name of section to generate
    :param existing_sections: List of sections to add to if a dependency or None if this is an outer-level abstract task
    """
    if not os.path.exists:
        existing_data = NEW_CONFIG_DATA
    else:
        existing_data = yaml.load(open(file_path, "r"), Loader=yaml.FullLoader)
        if existing_sections is not None:
            for section in existing_sections:
                assert section in existing_data.keys()
    # No sections to update - set as base level
    if existing_sections is None:
        existing_data[cfg_name] = CONFIG_OUTER_LEVEL
    # Else add to each section requested
    else:
        for section in existing_sections:
            existing_data[section]["dependencies"].update(CONFIG_INNER_LEVEL)
    yaml.dump(existing_data, open(os.path.splitext(file_path)[0] + ".tmp.yaml", "w"), Dumper=yaml.Dumper)


def create_class_file(file_path: str, class_name: str, cfg_name: str):
    """ Create class file at provided path. Will create given class using name and will
    register class as using cfg_name within config file

    :param file_path: Output path for class stub
    :param class_name: Name to give class
    :param cfg_name: Config section name
    """
    file_ptr = open(os.path.join(file_path, cfg_name.replace(".", "_") + ".py"), "w")
    file_ptr.write(BOILERPLATE.format(class_name, cfg_name, """{
                
            }"""))
    file_ptr.close()


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
            (("-s", "--config_sections"),
             {"help": "List all sections in existing config file to update. If not provided, will be added at "
                      "outer task level", "nargs": "+"})
        ),
        description="Create new EukMetaSanity class"
    )
    create_class_file(ap.args.path, ap.args.class_name, ap.args.config_name)
    if ap.args.config_path is not None:
        update_config_file(ap.args.config_path, ap.args.config_name, ap.args.config_sections)
