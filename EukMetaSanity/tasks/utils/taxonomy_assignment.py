"""
Module holds simple logic for passing taxonomy assignment/information between tasks
"""

from collections import namedtuple
from typing import Optional


class TaxonomyAssignment:
    """ This class holds the results of a taxonomic assignment by varying taxonomic levels

    """
    Assignment = namedtuple("Assignment", ("value", "tax_id", "score"))

    _tax_order = ["kingdom", "phylum", "_class", "order", "superfamily", "family", "genus", "species"]

    def __init__(self):
        """ Create empty assignment

        """
        self.kingdom: Optional[TaxonomyAssignment.Assignment] = None
        self.phylum: Optional[TaxonomyAssignment.Assignment] = None
        self._class: Optional[TaxonomyAssignment.Assignment] = None
        self.order: Optional[TaxonomyAssignment.Assignment] = None
        self.superfamily: Optional[TaxonomyAssignment.Assignment] = None
        self.family: Optional[TaxonomyAssignment.Assignment] = None
        self.genus: Optional[TaxonomyAssignment.Assignment] = None
        self.species: Optional[TaxonomyAssignment.Assignment] = None

    def assignment(self, level: str) -> Optional[Assignment]:
        """ Get Assignment object at given level.
        Returns None is level not found in file, or if level did not exist in file

        :param level: Taxonomy level string (e.g. kingdom, order, _class, etc.)
        :return: Assignment object or None
        """
        idx = TaxonomyAssignment._tax_order.index(level)
        _level = None
        while idx > 0 and _level is None:
            _level = getattr(self, TaxonomyAssignment._tax_order[idx], None)
            idx -= 1
        return _level

    def __repr__(self) -> str:
        """ Return string representation of assignment

        :return: Assignment as string
        """
        return "TaxonomyAssignment(kingdom={}, phylum={}, class={}, order={}, superfamily={}, family={}, genus={}, " \
               "species={})".format(
            *[
                str(v.value) if v is not None else "None" for v in (
                    self.kingdom, self.phylum, self._class, self.order, self.superfamily, self.family, self.genus,
                    self.species
                )
            ]
        )
