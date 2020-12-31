"""
Module contains PathManager that manages directory tree focused on using record_ids as subdirectories
"""

import os
from typing import List
from pathlib import Path
from plumbum import local


class PathManager:
    """ This class manages a set of working directories.
    Useful for keeping files organized
    :raises: AssertionError if base path does not exist
    """
    MAGS = "MAGS"

    def __init__(self, base_path: str):
        assert base_path is not None
        self._wdir = "wdir"
        self._base = str(Path(base_path).resolve())
        self._dbs = {}
        self._generate_directory_tree()

    # Get working dir
    @property
    def wdir(self):
        """ Working directory which PathManager object is maintaining.
        Multiple working directories are maintained within

        :return: str of full path to working directory
        """
        return os.path.join(self.base, self._wdir)

    # Get base dir
    @property
    def base(self):
        """ Base dir contains all working directories
        Typically it is the current directory, the directory the program was launched,
        or, possibly, some other user-specified directory (though this isn't yet implemented).

        :return: str of full path to base directory
        """
        return str(self._base)

    @property
    def dbs(self):
        """ Dictionary of subdirectories
        Subdirectories are maintained by a parent directory
        There are multiple parent directories within a given working directory

        :return: Dict of subdirectories and their paths
        """
        return self._dbs

    # Add record directory to wdir
    def add_dirs(self, record_id: str, _subdirs: List[str] = None):
        """ Add a subdirectory to a given record directory

        :param record_id: Directory within self.wdir to which to add subdirectories
        :param _subdirs: List of subdirectories to add to a directory within self.wdir
        :return:
        """
        assert record_id is not None
        # Record base dir
        local["mkdir"]["-p", os.path.join(self.wdir, str(record_id))]()
        # Additional dirs, if needed
        if _subdirs is not None:
            assert isinstance(_subdirs, list)
            for _subd in _subdirs:
                added_path = os.path.join(self.wdir, str(record_id), _subd)
                local["mkdir"]["-p", added_path]()
                self._dbs[_subd] = added_path

    # Get record directory in wdir
    def get_dir(self, record_id: str = None, subdir: str = None) -> str:
        """ Get full path to a subdirectory

        :param record_id: Directory in self.wdir with the label `record_id`
        :param subdir: Subdirectory within given label's directory
        :return: str of full path to record_id/subdir
        :raises: ValueError if unable to locate directory within storage
        """
        loc = self.wdir
        if record_id is not None:
            loc = os.path.join(loc, str(record_id))
        if subdir is not None:
            loc = self._dbs.get(subdir, "")
        if os.path.exists(loc):
            return loc
        raise ValueError(
            "%s%s dir does not exist in %s\n" % (record_id, ("".join(("/", subdir)) if subdir is not None else ""),
                                                 self.wdir)
        )

    # Generate initial directory tree
    def _generate_directory_tree(self):
        # Base output directory
        local["mkdir"]["-p", self._base]()
        # Working directory
        local["mkdir"]["-p", os.path.join(self._base, self._wdir)]()


if __name__ == "__main__":
    pass
