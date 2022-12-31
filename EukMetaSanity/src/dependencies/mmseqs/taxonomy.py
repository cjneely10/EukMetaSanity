import os
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import List, Union, Type, Dict, Optional, Tuple

from yapim import Task, DependencyInput, Result

TaxonomyResults = List[Tuple[str, Dict[str, Union[float, int, str]]]]


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
            result = {"taxid": node.data.ncbi_identifier,
                      "value": node.data.scientific_name,
                      "score": node.data.percent_mapped_reads}
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
                    stack.pop()
                    new_node = _Node.connect_nodes(stack.last().parent, node_info)
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
    def find_best_taxonomy(file: Path) -> TaxonomyResults:
        """
        Find the taxonomic label path to which the largest amount of reads mapped to the deepest taxonomic rank depths

        Calculates the path from root -> leaf to which the largest cumulative number of reads is mapped

        :param file: Result file from mmseqs taxonomyreport
        :return: Mapping consisting of {tax-level: {taxid: "", value: "", score: ""}}
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


class MMSeqsTaxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": os.path.join(self.wdir, f"{self.record_id}.txt"),
            "tax-db": Result(os.path.join(self.wdir, f"{self.record_id}-tax_db"))
        }
        tax_file = self.output["tax-report"]
        if os.path.exists(tax_file):
            self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(Path(tax_file).resolve())
        else:
            self.output["taxonomy"] = []

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        pass

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("MMSeqsCreateDB")]

    def run(self):
        # Search taxonomy db
        self.parallel(
            self.program[
                "taxonomy",
                self.input["MMSeqsCreateDB"]["db"],
                self.data[0],
                str(self.output["tax-db"]),
                os.path.join(self.wdir, "tmp"),
                (*self.added_flags),
                "--threads", self.threads,
                "--split-memory-limit", str(int(float(self.memory) * 0.7)) + "G",
            ]
        )
        # Generate taxonomy report
        self.single(
            self.program[
                "taxonomyreport",
                self.data[0],
                str(self.output["tax-db"]),
                self.output["tax-report"]
            ],
            "1:00:00"
        )
        self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(Path(self.output["tax-report"]))
