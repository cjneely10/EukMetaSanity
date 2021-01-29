"""
Module houses functionality to use in EukMS installation to download all required data dependencies for
official EukMS pipelines and run any merge/index generation steps needed for use in EukMS
"""
import os
from typing import Generator
from EukMetaSanity.api.data.data_types import Fasta, MMSeqsDB
from EukMetaSanity.api.data.mmseqs_operations import Merge, CreateTaxDB
from EukMetaSanity.api.data.mmseqs_index_types import CreateIndex, CreateLinIndex


def download_data(working_dir: str) -> Generator:
    """ Generate data download functors. Consumer should call each object in sequence.

    Current implementation:

    Download/generate mmseqs db for ODB

    Download/generate mmseqs db for MMETSP

    :param working_dir: Root directory for data downloads
    :return: Iterator over all official EukMS pipeline downloader functors
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    database_downloads = (
        Fasta(
            config_identifier="ortho_db",
            data_url="https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz",
            wdir=working_dir,
            expected="odb10v1_all_og_fasta.tab",
            unzip_command_args=["gunzip"],
        ),
        MMSeqsDB(
            config_identifier="mmetsp_db",
            data_url="https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/TAX_DBs/MMETSP/TaxDB_MMETSP.tar.gz",
            wdir=working_dir,
            expected="MMETSP",
            unzip_command_args=["tar", "-xzf", "-C", working_dir],
        ),
    )
    for database in database_downloads:
        yield database


def manage_downloaded_data(working_dir: str, create_index: bool, create_linindex: bool,
                           threads: int, split_mem_limit: str) -> Generator:
    """ Generate data utilities functors. Consumer should call each object in sequence.

    Current implementation:

    Generate ODB taxonomy database

    Merge ODB and MMETSP databases into single database

    :return: Iterator over all official EukMS download utility functors
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    fxns = [
        CreateTaxDB("mmseqs.input", working_dir, ["ortho_db"]),
        Merge("odb-mmetsp_db", ["ortho_db", "mmetsp_db"]),
    ]
    if create_index:
        fxns.append(CreateIndex(threads, split_mem_limit, ["ortho_db", "mmetsp_db", "odb-mmetsp_db"]))
    if create_linindex:
        fxns.append(CreateLinIndex(threads, split_mem_limit, ["ortho_db", "mmetsp_db", "odb-mmetsp_db"]))
    for fxn in fxns:
        yield fxn
