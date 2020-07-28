#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List, Dict
from collections import namedtuple
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse

GFFCoord = namedtuple("GFFCoord", ("evidence", "loc_type", "start", "end", "strand", "parent_id"))


class Coordinate:
    __slots__ = (
        "is_repeat_region",
        "is_exon",
        "loc_type",
        "parent_id",
        "evidence",
        "strand",
        "nucleotide",
    )

    def __init__(self, nucleotide: str):
        self.parent_id = None
        self.nucleotide = nucleotide
        self.is_exon = False
        if nucleotide == "N":
            self.is_repeat_region = True
        else:
            self.is_repeat_region = False
        self.evidence = set()
        self.strand = None
        self.loc_type = None

    def __repr__(self):
        # parent_id nucleotide is_exon is_intron is_repeat_region strand loc_type evidence
        return "%s %s %i %i %s %s %s" % (
            self.parent_id,
            self.nucleotide,
            self.is_exon,
            self.is_repeat_region,
            self.strand,
            self.loc_type,
            self.evidence
        )


# Convert gff3 file to dictionary
def gff3_to_dict(gff3_file: str) -> Dict[str, List[GFFCoord]]:
    gff3_fp = open(gff3_file, "r")
    out_data = defaultdict(list)
    for _line in gff3_fp:
        if _line.startswith("#"):
            continue
        line = _line.rstrip("\r\n").split()
        out_data[line[0]].append(
            GFFCoord(
                evidence=line[1],
                loc_type=line[2],
                start=int(line[3]),
                end=int(line[4]),
                strand=line[6],
                parent_id=line[8].split(";")[0].split("=")[1]
            )
        )
    return dict(out_data)


# Create list of Coordinate objects that describe each location on a chromosome
def generate_initial_region(record: SeqRecord) -> List[Coordinate]:
    return [
        Coordinate(val)
        for val in record.seq
    ]


def write_region(region: List[Coordinate], fp, record_id: str, _cds: List[SeqRecord], _min_seq_length: int, j: int):
    started = False
    start_pos = 0
    _id = ""
    end_pos = 0
    evidence = ""
    strand = ""
    for i, coord in enumerate(region):
        exon_list = []
        exon_started = False
        if coord.loc_type is not None and not started:
            started = True
            start_pos = i + 1
            evidence = ("metaeuk" if "metaeuk" in coord.evidence else "ab-initio")
            strand = coord.strand
        elif coord.loc_type is None and started:
            started = False
            end_pos = i + 1
            # Gather valid exon sections
            for k in range(start_pos - 1, end_pos):
                if region[k].is_exon and not exon_started:
                    exon_started = True
                    exon_list.append(k + 1)
                elif not region[k].is_exon and exon_started:
                    exon_started = False
                    exon_list[-1] = (exon_list[-1], k + 1)
        if end_pos - start_pos > _min_seq_length:
            gene_id = "gene" + str(j)
            # Transcript/gene info
            fp.write("".join((
                "\t".join((
                    record_id,
                    evidence,
                    "gene",
                    str(start_pos),
                    str(end_pos),
                    ".",
                    strand,
                    ".",
                    "ID=%s" % gene_id
                )),
                "\n",
            )))
            # Exon/CDS info
            for exon in exon_list:
                if not isinstance(exon, tuple):
                    continue
                fp.write("".join((
                    "".join((
                        "\t".join((
                            record_id,
                            evidence,
                            "exon",
                            str(exon[0]),
                            str(exon[1]),
                            ".",
                            strand,
                            "0",
                            "Parent=%s" % gene_id
                        )),
                        "\n"
                    )),
                )))
            # Add FASTA CDS if requested
            if len(exon_list) > 0:
                seq = Seq("".join((
                        char.nucleotide
                        for exon in exon_list for char in region[exon[0]:exon[1] + 1]
                    )))
            else:
                seq = Seq("".join((
                    char.nucleotide for char in region[start_pos - 1: end_pos]
                )))
            # Generate CDS sequence
            _cds.append(
                SeqRecord(
                    seq=seq,
                    id=record_id + "_" + gene_id,
                    description="strand=%s" % strand,
                    name="",
                )
            )
            # Reset counter/storage variables
            end_pos = 0
            j += 1
            evidence = ""
            strand = ""
    return j


# # Program driver logic

def _parse_args(ap: ArgParse):
    for _path in (ap.args.fasta_file, *ap.args.gff3_files):
        assert os.path.exists(_path)
    try:
        ap.args.min_seq_length = int(ap.args.min_seq_length)
    except ValueError as e:
        print(e)
        exit(1)
    assert ap.args.output_file is not None


def exonize(fasta_file: str, gff3_files: List[str], output_file: str, write_cds: str, write_prot: str, min_len: int):
    w = open(output_file, "w")
    # Get FASTA file as dict
    record_p = SeqIO.parse(fasta_file, "fasta")
    # Generate list of GFF data dictionaries
    gff_dict_list = [gff3_to_dict(_file) for _file in gff3_files]
    out_cds = []
    j = 1
    # Iterate over each record
    for record in record_p:
        # Create bare region
        region = generate_initial_region(record)
        for gff_dict in gff_dict_list:
            coord_datas = gff_dict.get(record.id, None)
            if coord_datas is None:
                continue
            for coord_data in coord_datas:
                # Parse region info
                for i in range(coord_data.start - 1, coord_data.end - 2):
                    region[i].parent_id = coord_data.parent_id
                    region[i].loc_type = coord_data.loc_type
                    region[i].evidence.add(coord_data.evidence)
                    region[i].strand = coord_data.strand
                    if not region[i].is_repeat_region:
                        region[i].is_exon = True
        # Write results in gff format
        j = write_region(region, w, record.id, out_cds, min_len, j)
    if write_cds is not None:
        SeqIO.write(out_cds, write_cds, "fasta")
    if write_prot is not None:
        SeqIO.write(find_orfs(out_cds), write_prot, "fasta")
    w.close()


# https://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence
def find_orfs(cds_list: List[SeqRecord]):
    for record in cds_list:
        longest = (0,)
        for nuc in (str(record.seq), str(record.reverse_complement().seq)):
            for m in re.finditer("ATG", nuc):
                pro = Seq(nuc)[m.start():].translate(to_stop=True)
                if len(pro) > longest[0]:
                    longest = (len(pro), m.start(), str(pro))
        if longest[0] >= 30:
            yield SeqRecord(
                seq=Seq(longest[2] + "*"),
                id=record.id,
                description=record.description
            )


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "Input FASTA file", "required": True}),
            (("-g", "--gff3_files"),
             {"help": "Input GFF3 files", "required": True, "nargs": "+"}),
            (("-o", "--output_file"),
             {"help": "Output path, default stdout"}),
            (("-c", "--cds"),
             {"help": "Output CDS sequences to path"}),
            (("-p", "--prot"),
             {"help": "Output amino acid sequences to path"}),
            (("-m", "--min_seq_length"),
             {"help": "Minimum CDS sequence to report with -c, default 500", "default": "500"}),
        ),
        description="Convert EukMetaSanity .merged.gff3 into exonized .nr.gff3 file"
    )
    _parse_args(_ap)
    exonize(
        _ap.args.fasta_file,
        _ap.args.gff3_files,
        _ap.args.output_file,
        _ap.args.cds,
        _ap.args.prot,
        _ap.args.min_seq_length
    )
