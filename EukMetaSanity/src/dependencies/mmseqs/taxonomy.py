import os
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import List, Union, Type, Dict, Optional

from yapim import Task, DependencyInput, Result


@dataclass
class _NodeInfo:
    """Taxonomy node representing a line in mmseqs taxonomy-report output file"""
    percent_mapped_reads: float
    number_reads_covered_by_clade: int
    number_reads_assigned_to_taxon: int
    rank: Optional[str]
    ncbi_identifier: int
    scientific_name: str

    @cached_property
    def scientific_name_cleaned(self) -> str:
        return self.scientific_name.lstrip(" ")

    @cached_property
    def calculate_level(self) -> int:
        for i, s in enumerate(self.scientific_name):
            if s != " ":
                return i // 2  # Each new level is formed by adding 2 spaces to the preceding level
        return 0


@dataclass
class _Node:
    """Parsed information from mmseqs/kraken-like output and connections to remainder of tree"""
    data: Optional[_NodeInfo]
    children: List["_Node"]
    parent: Optional["_Node"]

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def add_child(self, child: "_Node"):
        self.children.append(child)


class _NodeStack:
    def __init__(self):
        self._stack: List[_Node] = []

    def push(self, v: _Node):
        self._stack.append(v)

    def pop(self) -> Optional[_Node]:
        if self.is_empty():
            return None
        return self._stack.pop()

    def __len__(self) -> int:
        return len(self._stack)

    def is_empty(self) -> bool:
        return len(self) == 0

    def last(self) -> Optional[_Node]:
        if self.is_empty():
            return None
        return self._stack[-1]

    def clear(self):
        self._stack.clear()


class MMSeqsTaxonomyReportParser:
    _tax_order = ["superkingdom", "kingdom", "phylum", "class", "order", "superfamily", "family", "genus", "species"]

    @staticmethod
    def _parse_line(line: List[str]) -> _NodeInfo:
        """
        Create _NodeInfo from line that matches mmseqs/kraken specs

        :param line: File line split by tab
        :return: Parsed _NodeInfo
        :raises ValueError: If a portion of the line does not match mmseqs/kraken specifications
        """
        assert len(line) == 6, "file does not match mmseqs/kraken documentation specifications"
        return _NodeInfo(
            float(line[0]),
            int(line[1]),
            int(line[2]),
            None if line[3] == "no rank" else line[3],
            int(line[4]),
            line[5])

    @staticmethod
    def create_tree(file: Path) -> _Node:
        """
        Convert report file to taxonomy tree

        :param file: Result file from mmseqs taxonomy-report
        :raises ValueError: If parsing a field in the taxonomy file line fails (i.e., int(line[0]) fails, etc.)
        :return: Root node to taxonomy tree
        """
        stack = _NodeStack()
        root = _Node(None, [], None)
        stack.push(root)
        current_level: int = -1
        with open(file, "r") as file_ptr:
            for _line in file_ptr:
                line: List[str] = _line.rstrip("\r\n").split("\t")
                node_info: _NodeInfo = MMSeqsTaxonomyReportParser._parse_line(line)
                node_level: int = node_info.calculate_level
                # Node is a sibling of the current node
                if current_level == node_level:
                    parent = stack.last().parent
                    new_node = _Node(node_info, [], parent)
                    if parent is not None:
                        parent.add_child(new_node)
                # Node is a child of the current node
                elif node_level > current_level:
                    new_node = _Node(node_info, [], stack.last())
                    stack.last().add_child(new_node)
                    current_level = node_level
                # Node is attached to *an ancestor* of the current node
                else:
                    # Pop off from stack until `current_level` matches node_level
                    while current_level != node_level:
                        stack.pop()
                        current_level -= 1
                    new_node = _Node(node_info, [], stack.last().parent)
                stack.push(new_node)
        # Remove trailing references
        stack.clear()
        return root

    @staticmethod
    def find_best_taxonomy(file: Path, cutoff: float) -> Dict[str, Dict[str, Union[float, int, str, None]]]:
        """
        Find the taxonomic label path to which the largest amount of reads mapped to the deepest depths

        :param file: Result file from mmseqs taxonomy-report
        :param cutoff: Minimum % (0-100) of mapped reads that should be considered in generated tree
        :return: Mapping consisting of {tax-level: {taxid: "", value: "", score: ""}}
        :raises AssertionError: If line in taxonomy file is of improper length, or if function parameters are invalid
        """
        assert file.exists()
        assert 0 < cutoff < 100
        out = {level: {"taxid": None, "value": None, "score": None} for level in MMSeqsTaxonomyReportParser._tax_order}
        try:
            root = MMSeqsTaxonomyReportParser.create_tree(file)
        except ValueError:
            return out

        return out


class MMSeqsTaxonomy(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "tax-report": os.path.join(self.wdir, f"{self.record_id}.txt"),
            "taxonomy": {},
            "tax-db": Result(os.path.join(self.wdir, f"{self.record_id}-tax_db"))
        }
        tax_file = self.output["tax-report"]
        if os.path.exists(tax_file):
            self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(
                Path(tax_file).resolve(), float(self.config["cutoff"]))

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
        self.output["taxonomy"] = MMSeqsTaxonomyReportParser.find_best_taxonomy(
            Path(self.output["tax-report"]).resolve(), float(self.config["cutoff"]))
