"""
Base classes in use
"""
# pylint: disable=too-few-public-methods
import os
from pathlib import Path
from abc import abstractmethod
from typing import List, Optional, Set, Sequence
from plumbum import local
from plumbum.machines import LocalCommand
from EukMetaSanity import touch


class Command:
    """
    Class holds simple print-and-run command-line functionality
    """
    @staticmethod
    def run(cmd: LocalCommand) -> Optional[object]:
        """ Print and run local command

        :param cmd: Local plumbum command
        """
        print(cmd)
        return cmd()


class Data(Command):
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
        # Download data
        if not os.path.exists(self._data):
            self.run(local["wget"][data_url, "-O", self._data])
            # Unzip downloaded data
            self._unzip_data(unzip_command_args)
            touch(self._data)

    def __call__(self, *args, **kwargs):
        """ Parent class implementation of functor

        :param args:
        :param kwargs:
        :return:
        """

    def _unzip_data(self, unzip_command_args: Optional[List[str]]):
        """ Run command to unzip data based on list of passed strings

        :param unzip_command_args: Command strings to use to unzip data
        """
        if unzip_command_args is not None:
            self.run(local[unzip_command_args[0]][(*unzip_command_args[1:]), self._data])

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


class DataUtil(Command):
    """
    Parent class to all mmseqs data utilities. Defines functor that will handle specific utility operation.

    databases: set of database names under operation
    """
    def __init__(self, databases: Sequence[str]):
        """ Generate base class with member data

        :param databases: Databases under operation
        """
        self._databases: Set[str] = set(databases)

    @abstractmethod
    def __call__(self):
        """
        Child class implementation should perform given operation on stored database(s)
        """
        pass

    @property
    def databases(self) -> Set[str]:
        """ Get reference to set of stored databases

        :return: Set of stored database path/names
        """
        return self._databases
