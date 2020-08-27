#!/usr/bin/env python3
import os
import sqlite3
from collections import defaultdict

from Bio import SeqIO
from decimal import Decimal
from typing import List, Tuple, Generator
from EukMetaSanity.utils.arg_parse import ArgParse


def generate(db_path: str, columns: List[str]) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    # Create table
    conn.execute("""CREATE TABLE annotations (gene_id, %s)""" % ", ".join(columns))
    conn.commit()
    return conn


# Display summary
def summarize(conn: sqlite3.Connection, out_path: str, non_null: defaultdict):
    # Write results from database to ensure that it was written correctly
    cursor = conn.execute("SELECT * FROM annotations")
    fp = open(out_path, "w")
    # Write header
    keys = sorted(list(non_null.keys()))
    out_str = "# Total %s\n" % " ".join(("%s:%i" % (key, non_null[key]) for key in keys))
    fp.write(out_str)
    fp.write("".join((
        "\t".join([descr[0] for descr in cursor.description]),
        "\n"
    )))
    for row in cursor:
        fp.write("".join((
            "\t".join(row),
            "\n",
        )))
    fp.close()
    conn.close()
    print(out_str)


def _parse_args(_ap: ArgParse) -> List[Tuple[str, str]]:
    _file: str
    annotation_files: List[Tuple[str, str]] = []
    # Confirm valid input paths
    assert _ap.args.fasta_file is not None and os.path.exists(_ap.args.fasta_file)
    # Ensure number types are valid
    try:
        _ap.args.max_evalue = Decimal(_ap.args.max_evalue)
    except ValueError as e:
        print(e)
        exit(1)
    # Convert annotation string to tuple of (file_type, path)
    annot: str
    f_string: List[str]
    for annot in _ap.args.annotations:
        assert "=" in annot
        f_string = annot.split("=")
        assert os.path.exists(f_string[1]), f_string[1]
        annotation_files.append((f_string[0], f_string[1]))
    return annotation_files


def _kegg_iter(kegg_file: str) -> Generator[Tuple[str, str], str, None]:
    fp = open(kegg_file, "r")
    # Skip first 2 header lines
    all(next(fp) for _ in range(2))
    line: str
    _line: List
    for line in fp:
        # Only for KEGG verified hits
        if line[0] == "*":
            _line = line.rstrip("\r\n").split()
            yield _line[1], " ".join(_line[6:])


def _eggnog_iter(eggnog_file: str, max_evalue: Decimal) -> Generator[Tuple[str, str], Tuple[str, Decimal], None]:
    fp = open(eggnog_file, "r")
    # Skip first 4 header lines
    all(next(fp) for _ in range(4))
    for line in fp:
        if line[0] == "#":
            return
        _line = line.rstrip("\r\n").split()
        if Decimal(_line[2]) <= max_evalue:
            yield _line[0], _line[1]


def _mmseqs_iter(mmseqs_file: str, max_evalue: Decimal) -> Generator[Tuple[str, str], Tuple[str, Decimal], None]:
    fp = open(mmseqs_file, "r")
    for line in fp:
        _line = line.rstrip("\r\n").split()
        if Decimal(_line[10]) <= max_evalue:
            yield _line[0], _line[1]


def _fasta_id_iter(fasta_file: str) -> Generator[str, str, None]:
    data = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        data[str(record.id)] = len(record.seq)
    fp = open(fasta_file, "r")
    for line in fp:
        if line[0] == '>':
            _id = line[1:].rstrip("\r\n").split(" ")[0]
            if data[_id] >= 30:
                yield _id


def annotate(fasta_file: str, annotations: List[Tuple[str, str]], max_evalue: Decimal, out_prefix: str):
    # Create database, or load existing
    table_col_ids = sorted([annotation[0] for annotation in annotations])
    db_path = out_prefix + ".db"
    if os.path.exists(db_path):
        os.remove(db_path)
    conn = generate(
        db_path,
        table_col_ids
    )
    assert conn is not None
    # Insert blank row ids for values not already present
    rows = {}
    for fasta_id in _fasta_id_iter(fasta_file):
        rows[fasta_id] = [fasta_id, *["0" for _ in range(len(table_col_ids))]]
    _ids = set(rows.keys())
    i = 1
    non_null = 0
    # Build data based on passed files
    for annot_type, annot_path in annotations:
        if annot_type == "kegg":
            data_iter = _kegg_iter(annot_path)
        elif annot_type == "eggnog":
            data_iter = _eggnog_iter(annot_path, max_evalue)
        else:
            data_iter = _mmseqs_iter(annot_path, max_evalue)
        for record_id, record_annotation in data_iter:
            for val in ('"', "'"):
                record_annotation = record_annotation.replace(val, "")
            if record_id in _ids:
                rows[record_id][i] = record_annotation
        i += 1
    # Determine counts for number of annotations
    counts = defaultdict(int)
    for vals in rows.values():
        count = 0
        for i, val in enumerate(vals[1:]):
            if val != "0":
                count += 1
        for i in range(count, -1, -1):
            counts[i] += 1
    query_place_string = ",".join(("?" for _ in range(len(table_col_ids) + 1)))
    conn.executemany(
        "INSERT INTO annotations VALUES (%s)" % query_place_string, [(*val, ) for val in rows.values()]
    )
    conn.commit()
    summarize(conn, out_prefix + ".summary", counts)


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "FASTA file to annotate"}),
            (("-a", "--annotations"),
             {"help": "Annotation paths, in format type=file[ type2=file2]", "nargs": "+"}),
            (("-e", "--max_evalue"),
             {"help": "Max evalue, default 0.001", "default": "0.001"}),
            (("-s", "--summarize"),
             {"help": "Summarize db as tsv to path"}),
            (("-o", "--output"),
             {"help": "Output prefix, default out", "default": "out"}),
        ),
        description="Add annotations to FASTA/gff3 file from various results files"
    )
    _annotations = _parse_args(ap)
    annotate(ap.args.fasta_file, _annotations, ap.args.max_evalue, ap.args.output)
