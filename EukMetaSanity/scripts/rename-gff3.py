#!/usr/bin/env python3
import os
from typing import Dict
from EukMetaSanity.utils.arg_parse import ArgParse


def _parse_args(ap: ArgParse):
    for _arg in (ap.args.gff3_file, ap.args.id_mapping_file):
        assert os.path.exists(_arg)


def rename_gff3(gff3_file: str, id_mapping_file: str, output_file: str, use_id: str, reverse: bool):
    # Load ids into dict
    _file_prefix = (use_id if use_id is not None else os.path.basename(os.path.splitext(gff3_file)[0]))
    _ids_dict = _parse_ids_file(
        id_mapping_file,
        _file_prefix,
        reverse
    )
    # Write file
    out_fp = open(output_file, "w")
    in_fp = open(gff3_file, "r")
    for line in in_fp:
        # Write comments
        if line.startswith(_file_prefix[0]):
            _line = line.rstrip("\r\n").split("\t", maxsplit=1)
            out_fp.write(line.replace(_line[0], _ids_dict[_line[0]]))
        else:
            out_fp.write(line)

    out_fp.close()


def _parse_ids_file(ids_file: str, file_prefix: str, reverse: bool) -> Dict[str, str]:
    _ids_fp = open(ids_file, "r")
    line = next(_ids_fp)
    # Determine which id column to use
    if reverse:
        key, val = 0, 1
    else:
        key, val = 1, 0
    # Move to starting position in file
    while not line.startswith(file_prefix):
        line = next(_ids_fp)
    # Load into dictionary
    out = {}
    line = next(_ids_fp).rstrip("\r\n").split("\t")
    while True:
        if len(line) != 2:
            break
        out[line[key]] = line[val]
        try:
            line = next(_ids_fp).rstrip("\r\n").split("\t")
        except StopIteration:
            break
    return out


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("gff3_file",),
             {"help": "Output gff3 file from EukMetaSanity"}),
            (("id_mapping_file",),
             {"help": "Output ids.list file from EukMetaSanity"}),
            (("-o", "--output_file"),
             {"help": "Output new gff3 file to path, default stdout", "default": "/dev/stdout"}),
            (("-r", "--reverse"),
             {"help": "Insert EukMetaSanity ids into a file with its original ids",
              "default": False, "action": "store_true"}),
            (("-i", "--use_id"),
             {"help": "Use this id instead of basename of file"}),
        ),
        description="Rename .gff3 ids in EukMetaSanity output using ids from run"
    )
    _parse_args(_ap)
    rename_gff3(_ap.args.gff3_file, _ap.args.id_mapping_file, _ap.args.output_file, _ap.args.use_id, _ap.args.reverse)
