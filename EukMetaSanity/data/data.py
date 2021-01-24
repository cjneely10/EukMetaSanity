"""
Functionality to download and extract databases. Provide additional instructions for generating input files needed in
parsing functions like mmseqs createseqtaxdb
"""
import os
from pathlib import Path
from typing import Dict, List, Optional
from plumbum import local


class Data:
    """
    Base class for data downloaded. Provides functor to unzip data with provided command
    """

    def __init__(self, config_identifier: str, data_url: str, wdir: str, unzip_command_args: Optional[List[str]] = None):
        """ Base class constructor handles creating references to common member data

        Per API, this should mimic the following:

        data: /path/to/mmsetsp_db

        :param config_identifier: Name assigned in config file, like mmetsp_db
        :param data_url: url/ftp link for data download
        :param wdir: Working directory, if not provided defaults to os.getcwd()
        :param unzip_command_args: List of program name and arguments to run to unzip, ex. ["tar", "-xzf"]
        """
        self._config_identifier = config_identifier
        self._data_path = data_url
        if unzip_command_args is not None:
            local[unzip_command_args[0]]["", (*unzip_command_args[1:])]()
            self._unzipped = True
        else:
            self._unzipped = False
        self._wdir = Path(wdir).resolve()

    @property
    def wdir(self) -> Path:
        """ Working directory

        :return: Path to working directory
        """
        return self._wdir

    @property
    def config_name(self) -> str:
        """ Get database name assigned in config file.

        :return: Database name
        """
        return self._config_identifier

    @property
    def url(self) -> str:
        """ Get url of file downloaded to generate database

        :return: URL/FTP link of data downloaded
        """
        return self._data_path

    def data(self, expected: Optional[str] = None) -> Path:
        """ Get expected file path after unzipping

        :param expected: Provide expected file name, or will attempt to guess
        :raises: FileExistsError if unable to properly guess output data
        :return: Path to downloaded, possibly extracted, data
        """
        # Provided data, use expected name and working directory
        if expected is not None:
            output_path = os.path.join(self._wdir, expected)
        else:
            # Not provided - "guess" by removing extension of provided file and .tar
            output_path = os.path.join(self._wdir or os.getcwd(), os.path.basename(self._data_path))
            if self._unzipped:
                output_path = os.path.splitext(output_path)[0].replace(".tar", "")
        if not os.path.exists(output_path):
            raise FileNotFoundError(output_path)
        return Path(output_path).resolve()


class Fasta(Data):
    """
    Class represents a FASTA file data type download
    """

    def __init__(self, create_index: bool, create_linindex: bool, *args, **kwargs):
        """ FASTA format will be input into mmseqs createdb

        :param fasta_file: Path to FASTA file downloaded
        :param create_index: Set if a database index should be created
        :param create_linindex: Set if a database linindex should be created
        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)
        self._index = create_index
        self._linindex = create_linindex
        self._fasta_file = ""

    def __call__(self, *args, **kwargs):
        """
        Generate mmseqs database from a FASTA file
        """
        local["mmseqs"]["createdb", self.data(), self.config_name]()


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
