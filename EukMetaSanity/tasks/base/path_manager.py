"""
Module contains PathManager that manages directory tree focused on using record_ids as subdirectories
"""

import os
from typing import List, Dict
from pathlib import Path
from plumbum import local


class PathManager:
    """
    This class manages a set of working directories. Useful for keeping files organized
    """
    MAGS = "MAGS"

    def __init__(self, base_path: str):
        """ Create PathManager object rooted at `base_path`

        :param base_path: Base path to use to generate root of directory tree
        :raises: AssertionError if base path is not string or is null-string
        """
        assert isinstance(base_path, str) and len(base_path) > 0
        self._wdir = "wdir"
        self._base = str(Path(base_path).resolve())
        self._dbs = {}
        self._generate_directory_tree()

    @property
    def wdir(self) -> str:
        """ Working directory which PathManager object is maintaining.
        Multiple working directories are maintained within

        :return: str of full path to working directory
        """
        return os.path.join(self.base, self._wdir)

    @property
    def base(self) -> str:
        """ Base dir contains all working directories
        Typically it is the current directory, the directory the program was launched,
        or, possibly, some other user-specified directory (though this isn't yet implemented).

        :return: str of full path to base directory
        """
        return self._base

    @property
    def dbs(self) -> Dict[str, str]:
        """ Dictionary of subdirectories
        Subdirectories are maintained by a parent directory
        There are multiple parent directories within a given working directory

        :return: Dict of subdirectories and their paths
        """
        return self._dbs

    def add_dirs(self, record_id: str, _subdirs: List[str] = None):
        """ Add a subdirectory to a given record directory

        :param record_id: Directory within self.wdir to which to add subdirectories
        :param _subdirs: List of subdirectories to add to a directory within self.wdir
        :raises: AssertionError for improper types
        """
        assert isinstance(record_id, str) and len(record_id) > 0
        # Record base dir
        local["mkdir"]["-p", os.path.join(self.wdir, str(record_id))]()
        # Additional dirs, if needed
        if _subdirs is not None:
            assert isinstance(_subdirs, list)
            # Create all subdirectories provided in list and track them
            for _subd in _subdirs:
                added_path = os.path.join(self.wdir, str(record_id), _subd)
                local["mkdir"]["-p", added_path]()
                self._dbs[_subd] = added_path

    def get_dir(self, record_id: str = None, subdir: str = None) -> str:
        """ Get full path to a subdirectory. Returns working directory if no args passed

        :param record_id: Directory in self.wdir with the label `record_id`
        :param subdir: Subdirectory within given label's directory
        :raises: ValueError if unable to locate directory within storage
        :return: str of full path to record_id/subdir
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

    def _generate_directory_tree(self):
        """ Generate initial working directory tree within manager's root/base directory

        """
        # Base output directory
        local["mkdir"]["-p", self._base]()
        # Working directory
        local["mkdir"]["-p", os.path.join(self._base, self._wdir)]()
