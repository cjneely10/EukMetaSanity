import os
from numba import types
from pathlib import Path
from plumbum import local
from numba.experimental import jitclass

mkdir = local["mkdir"]


@jitclass((
        ("_base", types.string),
        ("_wdir", types.string),
))
class PathManager:

    def __init__(self, base_path):
        self._wdir = "wdir"
        self._base = str(Path(base_path).resolve())
        self._generate_directory_tree()

    # Get working dir
    def wdir(self):
        return os.path.join(self._base, self._wdir)

    # Get base dir
    def base(self):
        return str(self._base)

    # Add record directory to wdir
    def add_dir(self, record_id, _subdirs=None):
        # Record base dir
        mkdir["-p", os.path.join(self.wdir(), str(record_id))]()
        # Additional dirs, if needed
        if _subdirs is not None:
            if not isinstance(_subdirs, list):
                _subdirs = [_subdirs]
            all(mkdir["-p", os.path.join(self.wdir(), str(record_id), _subd)]() for _subd in _subdirs)

    # Get record directory in wdir
    def get_record_dir(self, record_id):
        loc = os.path.join(self.wdir(), str(record_id))
        if os.path.exists(loc):
            return loc
        raise ValueError("%s dir does not exist\n" % record_id)

    # Generate initial directory tree
    def _generate_directory_tree(self):
        # Base output directory
        mkdir["-p", self._base]()
        # Working directory
        mkdir["-p", os.path.join(self._base, self._wdir)]()
