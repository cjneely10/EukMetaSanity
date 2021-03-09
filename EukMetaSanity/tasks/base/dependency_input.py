"""
Module holds DependencyInput `struct` that is used to in TaskList.depends member
"""

from typing import Optional, Tuple


# pylint: disable=too-few-public-methods
class DependencyInput:
    """ Class handles parsing input for dependency requirements

    """
    def __init__(self, name: str, primary_input: str = "root", id_mapping=None):
        """ Create object

        :param name: Name of dependency, example mmseqs.createdb
        :param primary_input: Input to use, default is root.fna
        :param id_mapping: Specific expected_key:available_key to override within implementation
        """
        self._name = name
        self._input = primary_input
        self._id_mapping = None
        if id_mapping is not None:
            self._id_mapping = tuple(id_mapping)

    @property
    def name(self) -> str:
        """ Get name of dependency

        :return: str of name
        """
        return self._name

    @property
    def input(self) -> str:
        """ Get dependency input object string

        :return: Parsed input object string
        """
        return self._input

    @property
    def id_mapping(self) -> Optional[Tuple]:
        """ Get mapping to generic-level dependency input

        :return: Mapping provided at Task level
        """
        return self._id_mapping

    def __str__(self):
        return f"Dependency <Name: {self._name} Root: {self.input} ID-mapping: {self._id_mapping}>"
