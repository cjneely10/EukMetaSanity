#!/usr/bin/env python3

"""
Module holds logic to generate TaskList/Task class stub and associated config file sections
"""

import os
from EukMetaSanity.utils.arg_parse import ArgParse


BOILERPLATE = """import os
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import ProcessExecutionError, CommandNotFound
from EukMetaSanity import InvalidPathError, MissingDataError, InvalidProtocolError
from EukMetaSanity import program_catch, prefix, touch, set_complete


class {0}Iter(TaskList):
    name = "{1}"
    requires = []
    depends = []

    class {0}(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {2}

        @program_catch
        def run(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__({0}Iter.{0}, {0}Iter.name, *args, **kwargs)
"""

# TODO(1): Handle config file generation


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
        ),
        description="Create new EukMetaSanity class"
    )
    create_class_file(ap.args.path, ap.args.class_name, ap.args.config_name)
