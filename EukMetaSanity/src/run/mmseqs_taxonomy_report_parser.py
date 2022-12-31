from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import List, Tuple, Dict, Union, Optional

TaxonomyAssignment = Dict[str, Union[float, int, str]]
RankResult = Tuple[str, TaxonomyAssignment]
TaxonomyResults = List[RankResult]


def _create_taxonomy_result(taxid: int, scientific_name: str, score: float) -> TaxonomyAssignment:
    return {"taxid": taxid, "value": scientific_name, "score": score}


def _create_taxonomic_ranks() -> List[str]:
    _ranks = ["kingdom", "phylum", "class", "cohort", "order",
              "family", "tribe", "genus", "section", "species", "strain"]
    out = []
    for rank in _ranks[:-2]:
        out.append(f"super{rank}")
        out.append(rank)
        out.append(f"sub{rank}")
        out.append(f"infra{rank}")
    out.append(_ranks[-2])
    out.append(_ranks[-1])
    return out


# Ranks supported by current taxonomic databases


_supported_taxonomic_ranks = _create_taxonomic_ranks()


@dataclass
class _TaxonomyInfo:
    """Taxonomy node representing a line in mmseqs taxonomyreport output file"""
    percent_mapped_reads: float  # 0
    number_reads_assigned_to_taxon: int  # 2
    rank: Optional[str]  # 3
    ncbi_identifier: int  # 4
    scientific_name_string: str  # 5

    @cached_property
    def scientific_name(self) -> str:
        return self.scientific_name_string.lstrip(" ")

    @cached_property
    def calculate_level(self) -> int:
        for i, s in enumerate(self.scientific_name_string):
            if s != " ":
                return i // 2  # Each new level is formed by adding 2 spaces to the preceding level
        return 0


@dataclass
class _Node:
    """Parsed information from mmseqs/kraken-like output and connections to remainder of tree"""
    data: Optional[_TaxonomyInfo]
    children: List["_Node"]
    parent: Optional["_Node"]
    cumulative_read_count: int

    def add_child(self, child: "_Node"):
        if self.data is not None:
            child.cumulative_read_count += self.data.number_reads_assigned_to_taxon
        self.children.append(child)

    @staticmethod
    def connect_nodes(parent: Optional["_Node"], child_info: "_TaxonomyInfo") -> "_Node":
        child_node = _Node(child_info, [], parent, child_info.number_reads_assigned_to_taxon)
        if parent is not None:
            parent.add_child(child_node)
        return child_node

    def collect_taxonomy(self) -> TaxonomyResults:
        node = self
        out = []
        while node is not None and node.data is not None and node.data.scientific_name != "root":
            result = _create_taxonomy_result(node.data.ncbi_identifier,
                                             node.data.scientific_name,
                                             node.data.percent_mapped_reads)
            out.append((node.data.rank, result))
            node = node.parent
        out.reverse()
        return out


class _NodeStack:
    def __init__(self):
        self._stack: List[_Node] = []

    def push(self, v: _Node):
        self._stack.append(v)

    def pop(self) -> Optional[_Node]:
        if self.is_empty():
            return None
        return self._stack.pop()

    def last(self) -> Optional[_Node]:
        if self.is_empty():
            return None
        return self._stack[-1]

    def is_empty(self) -> bool:
        return len(self._stack) == 0

    def clear(self):
        self._stack.clear()


class MMSeqsTaxonomyReportParser:
    @staticmethod
    def _parse_line(line: List[str]) -> _TaxonomyInfo:
        """
        Create _NodeInfo from line that matches mmseqs/kraken specs

        :param line: File line split by tab
        :return: Parsed _NodeInfo
        :raises ValueError: If a portion of the line does not match mmseqs/kraken specifications and fails to parse
        """
        assert len(line) == 6, "file does not match mmseqs/kraken documentation specifications"
        return _TaxonomyInfo(float(line[0]), int(line[2]), line[3], int(line[4]), line[5])

    @staticmethod
    def _create_tree(file: Path) -> Tuple[_Node, List[_Node]]:
        """
        Convert report file to taxonomy tree

        :param file: Result file from mmseqs taxonomyreport
        :raises ValueError: If parsing a field in the taxonomy file line fails (i.e., int(line[0]) fails, etc.)
        :raises AssertionError: If file does not exist
        :return: Root node to taxonomy tree and a list of the tree's leaf nodes
        """
        assert file.exists()
        stack = _NodeStack()
        root = _Node(None, [], None, 0)
        stack.push(root)
        current_level: int = -1
        leaf_nodes: List[_Node] = []
        with open(file, "r") as file_ptr:
            for line in file_ptr:
                node_info: _TaxonomyInfo = MMSeqsTaxonomyReportParser._parse_line(line.rstrip("\r\n").split("\t"))
                node_level: int = node_info.calculate_level
                # Node is a sibling of the current node
                if current_level == node_level:
                    leaf_nodes.append(stack.last())
                    new_node: _Node = _Node.connect_nodes(stack.last().parent, node_info)
                    stack.pop()
                # Node is a child of the current node
                elif node_level > current_level:
                    current_level = node_level
                    new_node = _Node.connect_nodes(stack.last(), node_info)
                # Node is attached to *an ancestor* of the current node
                else:
                    leaf_nodes.append(stack.last())
                    # Pop off from stack until `current_level` matches node_level
                    while current_level != node_level:
                        stack.pop()
                        current_level -= 1
                    new_node = _Node.connect_nodes(stack.last().parent, node_info)
                    stack.pop()
                stack.push(new_node)
        # Add last node to leaf nodes
        last_stack_elem = stack.last()
        if last_stack_elem.data is not None:
            leaf_nodes.append(last_stack_elem)
        # Remove trailing references
        stack.clear()
        return root, leaf_nodes

    @staticmethod
    def _find_taxonomy(root: _Node, *taxonomic_ranks: str) -> bool:
        """
        Check if given taxonomic rank list corresponds to an existing path in the generated taxonomy tree

        :param root: Root of tree created by `_create_tree(...)`
        :param taxonomic_ranks: List of ranks to search as a path in the tree
        :return: Rank walk is present within the given tree
        """
        if len(taxonomic_ranks) == 0:
            return False
        for rank in taxonomic_ranks:
            rank_is_present = False
            for child in root.children:
                if child.data.scientific_name == rank:
                    root = child
                    rank_is_present = True
                    break
            if not rank_is_present:
                return False
        return True

    @staticmethod
    def _get_rank_results(results: TaxonomyResults, rank: str) -> Optional[RankResult]:
        """
        Find result for a given taxonomic rank

        :param results: Results from `cls.find_best_taxonomy`
        :param rank: Desired taxonomic rank
        :return: Result for rank, or None if rank was not present in results
        """
        # For clade, check for most-specific clade assignments
        if rank == "clade":
            _iter = results[::-1]
        else:
            _iter = results
        for rank_results in _iter:
            if rank_results[0] == rank:
                return rank_results
        return None

    @staticmethod
    def find_best_taxonomy(file: Path) -> TaxonomyResults:
        """
        Find the taxonomic label path to which the largest amount of reads mapped to the deepest taxonomic rank depths

        Calculates the path from root -> leaf to which the largest cumulative number of reads is mapped

        :param file: Result file from mmseqs taxonomyreport
        :return: List of results in format [(tax-rank, {taxid: "", value: "", score: ""})]
        :raises AssertionError: If line in taxonomy file is of improper length, or if function parameters are invalid
        :raises ValueError: If parsing a field in the taxonomy file line fails (i.e., int(line[0]) fails, etc.)
        """
        assert file.exists()
        root: _Node
        leaf_nodes: List[_Node]
        root, leaf_nodes = MMSeqsTaxonomyReportParser._create_tree(file)
        # Find leaf node with the highest cumulative read count
        leaf_nodes.sort(key=lambda node: node.cumulative_read_count, reverse=True)
        for leaf_node in leaf_nodes:
            if leaf_node.data is not None and leaf_node.data.scientific_name != "unclassified":
                return leaf_node.collect_taxonomy()
        # Failed to calculate taxonomy
        return []

    @staticmethod
    def find_assignment_nearest_request(taxonomy_results: TaxonomyResults, rank: str) -> RankResult:
        """
        Find assignment for a given taxonomic rank.

        If rank is not present, will attempt to search for closest taxonomic assignment upwards:

        For example, if rank="superfamily" and it is not present, this method may return the result for "infraorder",
        or any taxonomic rank that precedes superfamily.

        :param taxonomy_results: Results from `cls.find_best_taxonomy`
        :param rank: Taxonomic rank, e.g., "family"
        :return: Tuple consisting of rank (or nearest rank) and its result (tax-rank, {taxid: "", value: "", score: ""})
        """
        # First, check if exact rank is present
        matching_rank = MMSeqsTaxonomyReportParser._get_rank_results(taxonomy_results, rank)
        if matching_rank is not None:
            return matching_rank
        # If not present, find rank in supported taxonomic ranks and check for presence upwards
        try:
            rank_index = _supported_taxonomic_ranks.index(rank)
            for i in range(rank_index - 1, -1, -1):
                close_rank = MMSeqsTaxonomyReportParser._get_rank_results(taxonomy_results,
                                                                          _supported_taxonomic_ranks[i])
                if close_rank is not None:
                    return close_rank
        # Unsupported ranks will default to eukaryota assignment
        except ValueError:
            eukaryota_rank = MMSeqsTaxonomyReportParser._get_rank_results(taxonomy_results, "superkingdom")
            if eukaryota_rank is not None:
                return eukaryota_rank
        # If eukaryota was not assigned, will provide default assignment with `-1` mapping percent
        return "superkingdom", _create_taxonomy_result(2759, "Eukaryota", -1)
