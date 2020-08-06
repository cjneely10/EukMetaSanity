#!/usr/bin/env python3
import os
from typing import Dict
from EukMetaSanity.utils.arg_parse import ArgParse


def _parse_args(ap: ArgParse):
    for _arg in (ap.args.gff3_file, ap.args.id_mapping_file):
        assert os.path.exists(_arg)


def rename_gff3(gff3_file: str, id_mapping_file: str, output_file: str, use_id: str):

    _file_prefix = os.path.basename(os.path.splitext(gff3_file)[0])
    # Load ids into dict
    _ids_dict = _parse_ids_file(
        id_mapping_file,
        (use_id if use_id is not None else os.path.basename(os.path.splitext(gff3_file)[0]))
    )
    # Write file
    out_fp = open(output_file, "w")
    in_fp = open(gff3_file, "r")
    for line in in_fp:
        # Write comments
        if line.startswith("#"):
            out_fp.write(line)
        else:
            # Get first part of line
            line = line.split("\t", maxsplit=1)
            # Write replaced value
            out_fp.write("\t".join((
                _ids_dict[line[0]],
                line[1]
            )))
    out_fp.close()


def _parse_ids_file(ids_file: str, file_prefix: str, reverse: bool = False) -> Dict[str, str]:
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
    line = line.rstrip("\r\n").split("\t")
    while len(line) == 2:
        try:
            out[line[key]] = line[val]
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
    rename_gff3(_ap.args.gff3_file, _ap.args.id_mapping_file, _ap.args.output_file, _ap.args.use_id)
