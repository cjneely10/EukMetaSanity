"""
Functionality to download and extract databases. Provide additional instructions for generating input files needed in
parsing functions like mmseqs createseqtaxdb
"""
import os
from pathlib import Path
from typing import List, Optional
from plumbum import local


class Data:
    """
    Base class for data downloaded. Provides functor to unzip data with provided command
    """

    def __init__(self, config_identifier: str, data_url: str, wdir: str, expected: str,
                 unzip_command_args: Optional[List[str]] = None):
        """ Base class constructor handles creating references to common member data

        Per API, this should mimic the following:

        data: /path/to/mmsetsp_db

        :param config_identifier: Name assigned in config file, like mmetsp_db
        :param data_url: url/ftp link for data download
        :param wdir: Working directory, if not provided defaults to os.getcwd()
        :param expected: Download name for wget
        :param unzip_command_args: List of program name and arguments to run to unzip, ex. ["tar", "-xzf"]
        """
        self._config_identifier = config_identifier
        self._data_path = data_url
        self._wdir = Path(wdir).resolve()
        self._data = os.path.join(self._wdir, expected)
        self._unzip_data(unzip_command_args)
        # Download data
        local["wget", data_url, "-O", self._data]()

    def __call__(self, *args, **kwargs):
        """ Default implementation of functor does nothing

        :param args:
        :param kwargs:
        :return:
        """
        return

    def _unzip_data(self, unzip_command_args: Optional[List[str]]):
        """ Run command to unzip data based on list of passed strings

        :param unzip_command_args: Command strings to use to unzip data
        """
        if unzip_command_args is not None:
            local[unzip_command_args[0]]["", (*unzip_command_args[1:])]()
            self._unzipped = True
        else:
            self._unzipped = False

    @property
    def wdir(self) -> Path:
        """ Get working directory

        :return: Working directory
        """
        return self._wdir

    @property
    def url(self) -> str:
        """ Get url of file downloaded to generate database

        :return: URL/FTP link of data downloaded
        """
        return self._data_path

    @property
    def db_name(self) -> str:
        """ Get database name with working directory as joined path

        :return: Database name with working directory as joined path
        """
        return os.path.join(self._wdir, self._config_identifier)

    @property
    def data(self) -> Path:
        """ Get expected file path after unzipping

        :raises: FileExistsError if unable to properly guess output data
        :return: Path to downloaded, possibly extracted, data
        """
        # Provided data, use expected name and working directory
        if not os.path.exists(self._data):
            raise FileNotFoundError(self._data)
        return Path(self._data).resolve()


class Fasta(Data):
    """
    Class represents a FASTA file data type download
    """

    def __init__(self, *args, **kwargs):
        """ FASTA format will be input into mmseqs createdb

        :param fasta_file: Path to FASTA file downloaded
        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        """
        Generate mmseqs database from a FASTA file
        """
        local["mmseqs"]["createdb", self.data, self.db_name]()


class MMSeqsDB(Data):
    """
    Class represents an MMSEQs database download
    """

    def __init__(self, *args, **kwargs):
        """ Database will simply be extracted, default functor used

        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)


class MSA(Data):
    """
    Class represents a MSA to convert to MMseqs profile format
    """
    def __init__(self, *args, **kwargs):
        """ Format is FASTA and will convert to MMseqs profile

        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        """
        Overrides default call operator to create profile format
        """
        local["mmseqs"]["convertmsa", self.data, self.db_name + "-msa"]()
        local["mmseqs"]["msa2profile", self.db_name + "-msa", self.db_name]()

