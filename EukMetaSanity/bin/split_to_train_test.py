#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    USAGE = "split_to_train_test.py <input.fasta> <test-ratio=0.10>"
    assert len(sys.argv) == 3 or len(sys.argv) == 2, USAGE
    # Determine test ratio
    try:
        if len(sys.argv) == 3:
            sys.argv[2] = float(sys.argv[2])
        else:
            sys.argv.append(0.10)
    except ValueError as e:
        print(e)
        exit(1)


def split_train_test(input_fasta):
    test_out = []
    train_out = []
    record_fp = SeqIO.parse(input_fasta, "fasta")
    for num, record in enumerate(record_fp):
        print("%s/6" % str(num + 1))
        len_seq = len(record.seq)
        test_seq = "".join(val for i, val in enumerate(record.seq) if 0 <= i % 100000 <= 10000)
        train_seq = "".join(val for i, val in enumerate(record.seq) if 10000 < i % 100000 < 100000)
        # Store in lists to write
        for _list, _seq in zip((train_out, test_out), (train_seq, test_seq)):
            _list.append(
                SeqRecord(
                    id=record.id,
                    seq=Seq(_seq),
                    description=""
                )
            )
        for _text, _seq in zip(("Train", "Test"), (len(train_seq), len(test_seq))):
            print("%s len: %s" % (_text, str(_seq)))
        print("Total: %s, Total length: %s" % (len(train_seq) + len(test_seq), len_seq))
    # Write train/test results
    prefix = os.path.splitext(input_fasta)[0]
    for _out, _path in zip((train_out, test_out), (".train.fasta", ".test.fasta")):
        SeqIO.write(_out, prefix + _path, "fasta")


def main(argv):
    split_train_test(argv[1])


if __name__ == "__main__":
    parse_args()
    main(sys.argv)
