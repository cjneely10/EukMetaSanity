"""
Module contains helpers for dynamically loading modules
"""
import os
from inspect import isclass
from pkgutil import iter_modules
from pathlib import Path
from importlib import import_module


def get_modules(calling_file: str) -> dict:
    """ Dynamically load all modules contained at sublevel

    :param calling_file: __file__ from location where this function is called
    :return: Dict of task.name: Task object class
    """
    # iterate through the modules in the current package
    package_dir = Path(calling_file).resolve().parent
    out = {}
    for (_, module_name, _) in iter_modules([package_dir]):
        # import the module and iterate through its attributes
        module = import_module(
            "EukMetaSanity.tasks.{}.{}".format(
                os.path.splitext(os.path.basename(package_dir))[0],
                module_name
            )
        )
        for attribute_name in dir(module):
            attribute = getattr(module, attribute_name)
            if isclass(attribute):
                # Add the class to this package's variables
                out[attribute.name] = attribute
    return out
