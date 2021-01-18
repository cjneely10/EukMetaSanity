#!/usr/bin/env python3

# pylint: disable=invalid-name
"""
Script handles collecting results from all gene predictors and parsing into final output results
"""

import os
import re
import sys
from io import StringIO
from datetime import datetime
from operator import itemgetter
from collections import defaultdict
from typing import List, Dict, Tuple, Generator, Optional, Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from EukMetaSanity.utils.arg_parse import ArgParse

# Current identifiers that signify an abinitio-based prediction
ABINITIO_IDENTIFIERS = {"augustus", "GeneMark.hmm", "ab-initio"}


class Gene:
    """
    Class holds the data describing a given gene region. Merges ab-initio skeleton exons with evidence-based exons
    """
    # TODO: Parse in ab initio as List[List[exons]]
    def __init__(self, ab_initio_data: List, strand: str, term_exons: List, _tier: int):
        """ Create Gene using base ab-initio skeleton on a given strand. Define terminal exons in
        skeleton.

        :param ab_initio_data: List of initial exons
        :param strand: Strand of DNA
        :param term_exons: List of terminal exons
        """
        self.exons: List = ab_initio_data
        self.strand: str = strand
        self.num_ab_initio: int = len(ab_initio_data)
        # Number of ab-initio exons removed after confirming with evidence
        self.trimmed_ab_initio: int = len(ab_initio_data)
        # Number of additional exons found at evidence level
        self.added_evidence: int = 0
        # Set of terminal exons
        self.terminal_exons = set(term_exons)
        # Set tier quality of gene
        self._tier = _tier

    @property
    def tier(self) -> int:
        """ Get current tier of gene

        :return: Current tier stored, default is length of data initially passed
        """
        return self._tier

    @tier.setter
    def tier(self, t: int):
        """ Set current tier of gene

        :param t: Tier value
        """
        self._tier = t

    def filter(self, tier: int):
        """ Remove genes that do not have sufficient level of evidence. Set all gene-level metadata to 0 or empty
        collections

        :param tier: Minimum (inclusive) tier to keep in output
        :raises: AssertionError for 0 or negative-valued tier values
        """
        assert tier > 0
        if self._tier < tier:
            self.exons = []
            self.num_ab_initio = 0
            self.trimmed_ab_initio = 0
            self.added_evidence = 0
            self.terminal_exons = {}

    # pylint: disable=too-many-branches
    def tier0(self, evidence_data: List):
        """ Add list of exons from protein/transcriptomic-based evidence

        :param evidence_data: List of evidence-derived exons
        """
        if len(self.exons) == 0:
            return
        _len_ev, _len_ex = len(evidence_data), self.num_ab_initio
        if _len_ev >= _len_ex:
            _near = _len_ex / _len_ev
        else:
            _near = _len_ev / _len_ex
        # Remove ab initio exons if well-resolved protein mapping
        if _near >= .7:
            count = 1
            # Keep start exon
            out_exons = [self.exons[0]]
            end = len(self.exons)
            # Keep end if one exists
            if len(self.exons) > 1:
                out_exons.append(self.exons[-1])
                count += 1
                end -= 1
            # Search remaining for overlap
            for ab_exon in self.exons[1:end]:
                is_found = False
                _exon = None
                for exon in evidence_data:
                    # An overlap is found
                    if Gene.in_exon(exon, ab_exon):
                        is_found = True
                        _exon = exon
                        break
                if is_found:
                    count += 1
                    # Truncate exon to match evidence
                    if ab_exon in self.terminal_exons:
                        if self.strand == "+":
                            out_exons.append((ab_exon[0], _exon[1], ab_exon[2]))
                        else:
                            out_exons.append((_exon[0], ab_exon[1], ab_exon[2]))
                    else:
                        out_exons.append(ab_exon)
            self.trimmed_ab_initio = count
        else:
            out_exons = self.exons
        for exon in evidence_data:
            is_found = False
            for ab_exon in out_exons:
                if Gene.in_exon(ab_exon, exon):
                    is_found = True
                    break
            if not is_found:
                self.added_evidence += 1
                out_exons.append(exon)
        self.exons = out_exons
        if self.strand == "+":
            self.exons.sort(key=itemgetter(0))
        else:
            self.exons.sort(key=itemgetter(0), reverse=True)

    @staticmethod
    def in_exon(query_coord: Tuple[int, int, int], target_coord: Tuple[int, int, int]) -> bool:
        """ Returns if part of query coord overlaps target coord

        :param query_coord: Tuple for query location
        :param target_coord: Tuple for target location
        :return: Boolean if two regions overlap
        """
        return max(query_coord[0], target_coord[0]) <= min(query_coord[1], target_coord[1])


class GffReader:
    """
    Class reads a Gff file
    """
    def __init__(self, gff3_path: str):
        """ Open gff3 file

        :param gff3_path: Path to file
        :raises: AssertionError if gff3 file path does not exist
        """
        assert os.path.exists(gff3_path)
        self.fp = open(gff3_path, "r")
        self._count = 0

    def __iter__(self):
        """ Create iterator of genes in file

        :return: Iterator
        """
        return self.next_gene()

    # pylint: disable=stop-iteration-return
    def next_gene(self) -> Generator[defaultdict, None, None]:
        """ Get next gene in file

        :return: Generator over genes present in each section of gff3 file
        """
        _line: str
        line: List[str]
        line = next(self.fp).rstrip("\r\n").split("\t")
        while True:
            if line[0][0] == "#":
                line = next(self.fp).rstrip("\r\n").split("\t")
                continue
            self._count += 1
            # Putative gene
            gene_data = {
                "geneid": "%s_gene%i" % (line[0], self._count),
                "fasta-id": line[0],
                "strand": line[6],
            }
            transcripts = []
            line = next(self.fp).rstrip("\r\n").split("\t")
            while line[2] != "locus":
                # Read in transcript info - first line is a transcript: source,tstart,tend
                # Replace sources with `ab-initio`
                transcripts.append(
                    [(line[1].replace(val, "ab-initio") for val in ABINITIO_IDENTIFIERS), []]
                )
                line = next(self.fp).rstrip("\r\n").split("\t")
                # Add exon to current info
                while line[2] not in ("transcript", "locus", "gene"):
                    if line[2] == "CDS":
                        transcripts[-1][-1].append(
                            (int(line[3]), int(line[4]), int(line[7]))  # exstart,exend,offset
                        )
                    line = next(self.fp).rstrip("\r\n").split("\t")
            # Merge based on name
            data = defaultdict(list)
            terminal_exons = []
            for transcript in transcripts:
                data[transcript[0]].extend(transcript[-1])
                # Track the locations of terminal exons - e.g. these may have stop codons that will may need adjustment
                if transcript[0] in ABINITIO_IDENTIFIERS and len(transcript[-1]) > 1:
                    if gene_data["strand"] == "+":
                        terminal_exons.append(transcript[-1][-1])
                    else:
                        terminal_exons.append(transcript[-1][0])
            # Sort transcripts by position (and strand)
            for transcript in data.keys():
                if gene_data["strand"] == "+":
                    data[transcript].sort(key=itemgetter(0))
                else:
                    data[transcript].sort(key=itemgetter(0), reverse=True)
            # Store in data and yield
            gene_data["transcripts"] = data
            gene_data["terminal_exons"] = terminal_exons
            yield gene_data


class GffMerge:
    """
    Class loads all levels of evidence together and allows Gene class to merge together
    """
    def __init__(self, gff3_path: str, fasta_path: str, tier: int):
        """ Open gff3 file and fasta file. Load fasta data into memory

        :param gff3_path: Path
        :param fasta_path: Path
        :param tier: Parsing tier
        :raises: AssertionError if either path does not exist
        """
        assert os.path.exists(gff3_path) and os.path.exists(fasta_path)
        self.reader = GffReader(gff3_path)
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        self.tier = tier

    def merge(self) -> Generator[Tuple[Dict, SeqRecord, SeqRecord, List[int]], None, None]:
        """ Merge together genes from gff3 file with existing ab-initio skeleton

        :return: Iterator over merged genes
        """
        gene_data: Dict
        for gene_data in self.reader:
            # Generate initial exon structure
            gene = Gene(
                gene_data["transcripts"]["ab-initio"],
                gene_data["strand"],
                gene_data["terminal_exons"],
                len(gene_data["transcripts"])
            )
            # Tier 0 is conservative pairing of exons, removing exons without evidence and incorporating
            # exons that were not identified in ab initio predictions
            # TODO: Squash ab-initio exon tracks to single skeleton set
            if self.tier == 0:
                # Keep exons with evidence, add exons missed by ab-initio
                for val in gene_data["transcripts"].keys():
                    if val not in ABINITIO_IDENTIFIERS:
                        gene.tier0(gene_data["transcripts"][val])
                gene_data["transcripts"] = gene.exons
                yield (gene_data, *self.create_cds(gene_data, gene))
            else:
                # Add filter to genes that do not occur within user-defined threshold
                gene.filter(self.tier)
                # TODO: Calculate longest ORF if gene passes filter
                # TODO: At end, only return valid gene_data dict metadata with associated longest ORF CDS/AA
                gene_data["transcripts"] = gene.exons
                yield (gene_data, *self.create_cds(gene_data, gene))

    def create_cds(self, gene_data: dict, gene: Gene) -> Tuple[Optional[SeqRecord], Optional[SeqRecord], List[int]]:
        """ Create CDS from gene data in region

        :param gene_data: Gene metadata
        :param gene: Gene object to translate
        :return: CDS/protein records and list of offsets
        """
        if len(gene_data["transcripts"]) == 0:
            return None, None, []
        orig_seq = str(self.fasta_dict[gene_data["fasta-id"]].seq)
        strand = gene_data["strand"]
        out_cds: List[str] = []
        offsets: List[int] = []
        for exon in gene_data["transcripts"]:
            record = SeqRecord(seq=Seq(orig_seq[exon[0] - 1: exon[1]]))
            if strand == "-":
                record = record.reverse_complement()
            out_cds.append(str(record.seq))
            offsets.append(exon[2])
        if strand == "-":
            offsets.reverse()
        cds = Seq(GffMerge.longest_orf("".join(out_cds)))
        _prot_seq = cds.translate()
        _stats = "|".join(map(str, (
            gene.num_ab_initio,
            gene.trimmed_ab_initio,
            gene.added_evidence,
            len(gene.exons),
            len(_prot_seq),
        )))
        descr = "contig=%s strand=%s %s" % (
            gene_data["fasta-id"],
            strand,
            _stats
        )
        gene_data["stats"] = _stats
        return (
            SeqRecord(
                id=gene_data["geneid"],
                name="",
                description=descr,
                seq=_prot_seq
            ),
            SeqRecord(
                id=gene_data["geneid"],
                name="",
                description=descr,
                seq=cds
            ),
            offsets
        )

    @staticmethod
    def longest_orf(sequence: str) -> str:
        """ Find longest ORF in a sequence

        :param sequence: Str sequence to search
        :return: Longest ORF in sequence
        """
        longest = ""
        possible_starts = ("ATG", "CCA", "CTG", "CAG")
        for start in possible_starts:
            start_pos = (m.start() for m in re.finditer(start, sequence))
            for pos in start_pos:
                idx = set()
                orf = GffMerge.l_orf_helper(sequence, idx, StringIO(), 3, pos)
                if len(orf) > len(longest):
                    longest = orf
        return longest

    @staticmethod
    def l_orf_helper(sequence: str, idx: Set[int], out: StringIO, k: int, start: int) -> str:
        """ Recursive helper to search most possible combinations of codons in a string

        :param sequence: Sequence to search
        :param idx: Set of coordinates already searched
        :param out: CDS string being build from underlying region
        :param k: Offset to search
        :param start: Start position to search
        :return: Longest ORF
        """
        possible_ends = ("TAG", "TAA", "TGA")
        for offset in range(k):
            if start + offset + k <= len(sequence):
                if start + offset in idx:
                    continue
                pos = start + offset
                idx.add(pos)
                s = sequence[pos: pos + k]
                out.write(s)
                if s in possible_ends:
                    return out.getvalue()
                return GffMerge.l_orf_helper(sequence, idx, out, k, pos + k)
        return out.getvalue()


# pylint: disable=too-few-public-methods
class GffWriter:
    """
    Class writes results of merging all data and evidence
    """
    def __init__(self, in_gff3_path: str, fasta_file: str, output_prefix: str, tier: int):
        """ Create write object based on input data

        :param in_gff3_path: Initial gff3 file
        :param fasta_file: FASTA file associated with gff3 file
        :param output_prefix: Out prefix to write
        :param tier: Tiered output to pass to GffMerge object
        """
        self.in_fp = open(in_gff3_path, "r")
        self.base = output_prefix
        self.out_fp = open(self.base + ".nr.gff3", "w")
        self.merger = GffMerge(in_gff3_path, fasta_file, tier)

    def write(self):
        """ Write parsed results in gff3 format

        """
        self.out_fp.write("# EukMetaSanity annotations generated %s\n" % datetime.now().strftime("%Y/%m/%d %H:%M"))
        out_prots: List[SeqRecord] = []
        out_cds: List[SeqRecord] = []
        current_id = ""
        for gene_dict, prot, cds, offsets in self.merger.merge():
            if gene_dict["fasta-id"] != current_id:
                current_id = gene_dict["fasta-id"]
                self.out_fp.write("# Region %s\n" % current_id)
            gene = GffWriter._gene_dict_to_string(gene_dict, offsets)
            if gene is not None:
                self.out_fp.write(gene)
            if prot and len(prot.seq) > 0:
                out_prots.append(prot)
            if cds and len(cds.seq) > 0:
                out_cds.append(cds)
        SeqIO.write(out_prots, self.base + ".faa", "fasta")
        SeqIO.write(out_cds, self.base + ".cds.fna", "fasta")

    # pylint: disable=inconsistent-return-statements
    @staticmethod
    def _gene_dict_to_string(gene_data: Dict, offsets: List[int]) -> Optional[str]:
        """ Convert data in gene dict to string format

        :param gene_data: Collected gene data
        :param offsets: Offsets associated with each CDS region
        :return: String of gene data, if any is present
        """
        gene_data["transcripts"].sort(key=itemgetter(0))
        gene_id = gene_data["geneid"]
        mrna_id = gene_id + "-mRNA"
        version = "EukMS"
        ss = StringIO()
        if len(gene_data["transcripts"]) == 0:
            return
        ss.write("".join((
            "\t".join((
                gene_data["fasta-id"], version,
                "gene", str(gene_data["transcripts"][0][0]), str(gene_data["transcripts"][-1][1]),
                ".", gene_data["strand"], ".", "ID=%s;Scores=%s" % (gene_id, gene_data["stats"])
            )),
            "\n",
            "\t".join((
                gene_data["fasta-id"], version,
                "mRNA", str(gene_data["transcripts"][0][0]), str(gene_data["transcripts"][-1][1]),
                ".", gene_data["strand"], ".", "ID=%s;Parent=%s" % (mrna_id, gene_id)
            )),
            "\n",
        )))
        for j, exon_tuple in enumerate(gene_data["transcripts"], start=1):
            ss.write("".join((
                "\t".join((
                    gene_data["fasta-id"], version,
                    "exon", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], ".", "ID=%s-exon%i;Parent=%s" % (mrna_id, j, mrna_id)
                )),
                "\n",
                "\t".join((
                    gene_data["fasta-id"], version,
                    "CDS", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], str(offsets[j - 1]), "ID=%s-cds%i;Parent=%s" % (mrna_id, j, mrna_id)
                )),
                "\n",
            )))
        return ss.getvalue()


def parse_args(ap: ArgParse):
    """ Confirm input arguments are valid types and within valid ranges

    :param ap: ArgParse object
    :raises: AssertionError for improperly formatted data
    """
    try:
        ap.args.tier = int(ap.args.tier)
    except ValueError:
        print("Tier must be integer")
        raise AssertionError
    assert ap.args.tier >= 0, "Tier must be positive"

    try:
        ap.args.recursion_limit = int(ap.args.recursion_limit)
    except ValueError:
        print("Recursion limit must be integer")
        raise AssertionError
    assert ap.args.recursion_limit >= 0, "Recursion limit must be positive"

    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file), _file
    if ap.args.output_prefix is None:
        ap.args.output_prefix = os.path.splitext(ap.args.gff3_file)[0]


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("-g", "--gff3_file"),
             {"help": ".all.gff3 file", "required": True}),
            (("-f", "--fasta_file"),
             {"help": "FASTA file", "required": True}),
            (("-o", "--output_prefix"),
             {"help": "Output prefix, default is path/prefix of gff3_file"}),
            (("-t", "--tier"),
             {"help": "Tiered output, any value greater than 1, default is 0 for owned parsing", "default": 0}),
            (("-r", "--recursion_limit"),
             {"help": "Override python's default recursion limit, default 2000000", "default": "2000000"}),
        ),
        description="GFF3 output final annotations as <prefix>.nr.gff3"
    )

    parse_args(_ap)
    sys.setrecursionlimit(_ap.args.recursion_limit)
    writer = GffWriter(_ap.args.gff3_file, _ap.args.fasta_file, _ap.args.output_prefix, _ap.args.tier)
    writer.write()
