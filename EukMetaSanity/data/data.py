"""
Module contains logic to automate data download process on EukMetaSanity installation
"""
from pathlib import Path
from typing import Dict, List, Optional
from collections import namedtuple
from plumbum import local

UrlInfo = namedtuple("UrlInfo", ("url", "tar", "flags", "gz", "type"))


class Data:
    def __init__(self, db_name: str, data_path: str, unzip_command_args: Optional[List[str]] = None):
        self._db_name = db_name
        self._data_path = data_path
        self._unzip_command = unzip_command_args

    def __call__(self, *args, **kwargs):
        if self._unzip_command is not None:
            local[self._unzip_command[0]]["", (*self._unzip_command[1:])]()


class Fasta(Data):
    def __init__(self, fasta_file: Path, create_index: bool, create_linindex: bool, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._index = create_index
        self._linindex = create_linindex
        self._fasta_file = fasta_file

    def __call__(self, *args, **kwargs):
        super().__call__(*args, **kwargs)


def data_urls() -> Dict[str, UrlInfo]:
    """ Return dictionary of config_replacement_string: instructions for unpackaging data

    Upon first installation, all config files in the EukMetaSanity/config directory will have
    values corresponding to keys of the dictionary returned by this function replaced with their
    respective paths post-installation.

    This allows seamless linking of downloading required data for new pipelines in new
    installations

    :return: dict of str: UrlInfo consisting of replacement_string: instructions to follow to unpackage data
    """
    return {
        "ortho_db": UrlInfo(
            url="https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz",
            flags="",
            tar=False,
            gz=True,
            type="FASTA",
        ),
        "rfam_db": UrlInfo(
            url="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz",
            flags="",
            tar=False,
            gz=False,
            type="profile",
        ),
        "mmetsp_db": UrlInfo(
            url="https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/TAX_DBs/MMETSP/TaxDB_MMETSP.tar.gz",
            flags="-xzf",
            tar=True,
            gz=False,
            type="FASTA",
        ),
    }
