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


def main(argv):
	out_records = []
	w = open("/dev/stdout", "w")
	for _file in glob.glob(argv[1]):
		seq = ""
		record_id = os.path.basename(os.path.splitext(_file)[0])
		record_p = SeqIO.parse(_file, "fasta")
		for record in record_p:
			seq += str(record.seq).replace("N", "").replace("n", "").upper()
		out_records.append(
			SeqRecord(
				seq=Seq(seq),
				id=record_id,
				description="",
			)
		)
		w.write("".join([record_id, "\n"]))
	SeqIO.write(out_records, argv[2], "fasta")
	w.close()


if __name__ == "__main__":
	parse_args(sys.argv)
	main(sys.argv)
