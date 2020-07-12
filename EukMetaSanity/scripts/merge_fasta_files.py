#!/usr/bin/env python3
import os
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args(argv):
	USAGE = "Usage: merge_fasta_files.py <fasta-glob-statement> <output-filename>"
	# Check for help
	if "-h" in argv or "--help" in argv:
		print(USAGE)
		exit(0)
	# Confirm number
	assert len(argv) == 3, USAGE
	# Confirm existence of FASTA files
	for _file in glob.glob(argv[1]):
		assert os.path.exists(_file)


def write_merged(_files, output_path):
	SeqIO.write(
		merge(_files),
		output_path,
		"fasta"
	)


def merge(glob_statement):
	for _file in glob.glob(glob_statement):
		seq = ""
		record_id = os.path.basename(os.path.splitext(_file)[0])
		record_p = SeqIO.parse(_file, "fasta")
		for record in record_p:
			seq += str(record.seq).replace("N", "").replace("n", "").upper()
		yield SeqRecord(
			seq=Seq(seq),
			id=record_id,
			description="",
		)


def main(argv):
	write_merged(argv[1], argv[2])


if __name__ == "__main__":
	parse_args(sys.argv)
	main(sys.argv)
