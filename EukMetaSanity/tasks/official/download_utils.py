"""
Module houses functionality to use in EukMS installation to download all required data dependencies for
official EukMS pipelines and run any merge/index generation steps needed for use in EukMS
"""
import os
from typing import Generator
from EukMetaSanity.api.data.data_types import Fasta, MMSeqsDB
from EukMetaSanity.api.data.mmseqs_index_types import CreateIndex, CreateLinIndex
from EukMetaSanity.tasks.official.download_parsing_functions import odb_tax_parse
from EukMetaSanity.api.data.mmseqs_operations import ConcatDBs, CreateTaxDBs, CreateMappingFiles


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
            expected="odb10v1_all_og_fasta.tab.gz",
            unzip_command_args=["gunzip"],
        ),
        MMSeqsDB(
            config_identifier="mmetsp_db",
            data_url="https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/TAX_DBs/MMETSP/TaxDB_MMETSP.tar.gz",
            wdir=working_dir,
            expected="TaxDB_MMETSP.tar.gz",
            unzip_command_args=["tar", "-xvzf", "-C", os.path.basename(working_dir)],
        ),
    )
    for database in database_downloads:
        yield database


def parsing_operations(working_dir: str) -> Generator:
    """ Parse downloaded data, if needed, to create any mapping files that are used to build taxonomy databases

    :param working_dir: Directory containing downloaded data
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    parsing_fxns = (
        CreateMappingFiles(
            [odb_tax_parse],
            working_dir,
            [os.path.join(working_dir, "ortho_db")]
        ),
    )
    for parsing_fxn in parsing_fxns:
        yield parsing_fxn


def manage_downloaded_data(working_dir: str, create_index: bool, create_linindex: bool,
                           threads: int) -> Generator:
    """ Generate data utilities functors. Consumer should call each object in sequence.

    Current implementation:

    Generate ODB taxonomy database

    Merge ODB and MMETSP databases into single database

    :return: Iterator over all official EukMS download utility functors
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    fxns = [
        CreateTaxDBs(working_dir, ["ortho_db"]),
        ConcatDBs(working_dir, "odb-mmetsp_db", ["ortho_db", "MMETSP"]),
        CreateTaxDBs(working_dir, ["odb-mmetsp_db"]),
    ]
    if create_index:
        fxns.append(CreateIndex(threads, working_dir, ["ortho_db", "MMETSP", "odb-mmetsp_db"]))
    if create_linindex:
        fxns.append(CreateLinIndex(threads, working_dir, ["ortho_db", "MMETSP", "odb-mmetsp_db"]))
    for fxn in fxns:
        yield fxn
