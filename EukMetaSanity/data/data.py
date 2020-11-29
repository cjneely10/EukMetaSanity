from typing import Dict
from collections import namedtuple

"""
Class to handle tracking all required datasets for each task

"""

UrlInfo = namedtuple("UrlInfo", ("url", "tar", "flags", "gz", "type"))


def data_urls() -> Dict[str, UrlInfo]:
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


if __name__ == "__main__":
    pass
