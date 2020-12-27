from typing import Dict


class DependencyInput:
    """ Class handles parsing input for dependency requirements

    """
    def __init__(self, name: str, additional_data: Dict[str, object] = {}, primary_input: str = "root"):
        """ Create object

        :param name: Name of dependency, example mmseqs.createdb
        :param additional_data: Dict of additional data to pass to dependency
        :param primary_input: Input to use, default is root.fna
        """
        self._name = name
        self._args = additional_data
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

    @property
    def args(self) -> Dict[str, object]:
        """ Get additional passed data to dependency

        :return:
        """
        return self._args
