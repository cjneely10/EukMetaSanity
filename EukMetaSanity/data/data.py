"""
Module contains logic to automate data download process on EukMetaSanity installation
"""

from typing import Dict
from collections import namedtuple


# pylint: disable=pointless-string-statement
""" NamedTuple consisting of options for downloading/unpackaging data
Logic is present within download-data.py script to handle these specific types:
tar, gz, and the flags passed to them
FASTA and mmseqs profile databases
"""
UrlInfo = namedtuple("UrlInfo", ("url", "tar", "flags", "gz", "type"))


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
            type="profile"
        )
    }
