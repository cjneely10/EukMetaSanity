import os
from typing import List
from pathlib import Path
from plumbum import local

mkdir = local["mkdir"]


class PathManager:
    def __init__(self, base_path: str):
        assert base_path is not None
        self._wdir = "wdir"
        self._base = str(Path(base_path).resolve())
        self._dbs = {}
        self._generate_directory_tree()

    # Get working dir
    @property
    def wdir(self):
        return os.path.join(self.base, self._wdir)

    # Get base dir
    @property
    def base(self):
        return str(self._base)

    @property
    def dbs(self):
        return self._dbs

    # Add record directory to wdir
    def add_dirs(self, record_id: str, _subdirs: List[str] = None):
        assert record_id is not None
        # Record base dir
        mkdir["-p", os.path.join(self.wdir, str(record_id))]()
        # Additional dirs, if needed
        if _subdirs is not None:
            assert isinstance(_subdirs, list)
            for _subd in _subdirs:
                added_path = os.path.join(self.wdir, str(record_id), _subd)
                mkdir["-p", added_path]()
                self._dbs[_subd] = added_path

    # Get record directory in wdir
    def get_dir(self, record_id: str = None, subdir: str = None):
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
        mkdir["-p", self._base]()
        # Working directory
        mkdir["-p", os.path.join(self._base, self._wdir)]()


if __name__ == "__main__":
    pass
