#!/usr/bin/env python3
import os
import sqlite3
import string
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
def summarize(conn: sqlite3.Connection, out_path: str):
    fp = open(out_path, "w")
    # Write header
    cursor = conn.execute("SELECT * FROM annotations")
    fp.write("".join((
        "\t".join([descr[0] for descr in cursor.description]),
        "\n"
    )))
    found = 0
    for row in cursor:
        for _r in row[1:]:
            if _r != "0":
                found += 1
                break
        fp.write("".join((
            "\t".join(row),
            "\n",
        )))
    fp.close()
    conn.close()
    print("Found: %i" % found)


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
        assert os.path.exists(f_string[1])
        annotation_files.append((f_string[0], f_string[1]))
    return annotation_files


def _kegg_iter(kegg_file: str) -> Generator[Tuple[str, str], str, None]:
    fp = open(kegg_file, "r")
    # Skip first 2 header lines
    all(next(fp) for _ in range(2))
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
    fp = open(fasta_file, "r")
    for line in fp:
        if line[0] == '>':
            yield line[1:].rstrip("\r\n").split(" ")[0]


def annotate(fasta_file: str, annotations: List[Tuple[str, str]], max_evalue: Decimal):
    # Create database, or load existing
    table_col_ids = sorted([annotation[0] for annotation in annotations])
    db_path = "annotations.db"
    if os.path.exists(db_path):
        os.remove(db_path)
    conn = generate(
        db_path,
        table_col_ids
    )
    assert conn is not None
    # Insert blank row ids for values not already present
    rows = []
    for fasta_id in _fasta_id_iter(fasta_file):
        rows.append(
            (
                fasta_id,
                *["0" for _ in range(len(table_col_ids))],
            )
        )
    query_place_string = ",".join(("?" for _ in range(len(rows[0]))))
    conn.executemany("INSERT INTO annotations VALUES (%s)" % query_place_string, rows)
    conn.commit()
    # Add data
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
            conn.execute(
                "UPDATE annotations SET %s = '%s' WHERE gene_id = '%s'" % (
                    annot_type, record_annotation, record_id
                )
            )
        conn.commit()
    summarize(conn, os.path.splitext(db_path)[0] + ".tsv")


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "FASTA file to annotate"}),
            (("-a", "--annotations"),
             {"help": "Annotation paths, in format type=file[ type2=file2]", "nargs": "+"}),
            (("-e", "--max_evalue"),
             {"help": "Max evalue, default 1E-5", "default": "1E-5"}),
            (("-s", "--summarize"),
             {"help": "Summarize db as tsv to path"}),
        ),
        description="Add annotations to FASTA/gff3 file from various results files"
    )
    # Summarize
    if ap.args.summarize is not None:
        assert os.path.exists(ap.args.summarize)
        summarize(sqlite3.connect(ap.args.summarize), os.path.splitext(ap.args.summarize)[0] + ".tsv")
        exit(0)

    # Or build
    _annotations = _parse_args(ap)
    annotate(ap.args.fasta_file, _annotations, ap.args.max_evalue)
