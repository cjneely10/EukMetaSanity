from typing import Tuple


class DependencyInput:
    """ Class handles parsing input for dependency requirements

    """
    def __init__(self, name: str, primary_input: str = "root"):
        """ Create object

        :param name: Name of dependency, example mmseqs.createdb
        :param primary_input: Input to use, default is root.fna
        """
        self._name = name
        self._input = primary_input

    @property
    def name(self) -> str:
        """ Get name of dependency

        :return: str of name
        """
        return self._name

    @property
    def input(self) -> str:
        """ Get dependency input object string

        :return: Parsed input object string as tuple
        """
        return self._input
