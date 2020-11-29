import os
from typing import Dict

"""
General functions

"""


def touch(_path: str):
    """ Mimic touch command from linux
    """
    open(_path, "w").close()


def prefix(_path: str) -> str:
    """ For path /path/to/file.txt return file
    """
    return os.path.basename(os.path.splitext(_path)[0])
