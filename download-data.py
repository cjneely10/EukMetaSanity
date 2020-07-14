#!/usr/bin/env python3
import os
from plumbum import local
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.tasks.manager.data import data_urls
from EukMetaSanity.utils.path_manager import PathManager

# Dependencies
cp = local["cp"]
sed = local["sed"]
tar = local["tar"]
wget = local["wget"]
gunzip = local["gunzip"]
mmseqs = local["mmseqs"]


def _print_and_run(cmd):
    print(cmd)
    cmd()


def _parse_args(ap: ArgParse):
    if ap.args.index:
        assert ap.args.build, "Build must be set in order to index!"
    try:
        ap.args.threads = int(ap.args.threads)
    except ValueError as e:
        print(e)
        exit(1)
    return ap


# Download all data with wget from the data.py script - Part of API
def run(ap: ArgParse, out_dir: str):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Download each URL to folder
    for _id, url in data_urls().items():
        # Download
        new_dir = os.path.join(out_dir, _id)
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        _file = os.path.join(new_dir, os.path.basename(url.url))
        if not os.path.exists(_file) or ap.args.rewrite:
            _print_and_run(wget[url.url, "-O", _file])
            # Tar/gunzip
            if url.tar:
                _print_and_run(tar[url.flags, _file])
                _file = os.path.splitext(_file.replace(".tar", ""))[0]
            elif url.gz:
                _print_and_run(gunzip[_file])
                _file = os.path.splitext(_file)[0]
            # Generate MMseqs2 database
            _out = os.path.splitext(_file)[0] + "_db"
            if ap.args.build:
                # Generate database
                _print_and_run(
                    mmseqs[
                        "createdb",
                        _file,
                        _out
                    ]
                )
                # Add taxonomy info from ncbi
                _create_tax_db(_out, out_dir)
                # Generate linear index
                _print_and_run(
                    mmseqs[
                        "createlinindex",
                        _out,
                        os.path.join(os.path.basename(_out), "tmp"),
                        "--threads", str(ap.args.threads),
                        "--split-memory-limit", ap.args.max_mem,
                    ]
                )
            # Create MMseqs2 index file
            if ap.args.index:
                _print_and_run(
                    mmseqs[
                        "createindex",
                        _out,
                        os.path.join(os.path.basename(_out), "tmp"),
                        "--threads", str(ap.args.threads),
                        "--split-memory-limit", ap.args.max_mem,
                    ]
                )
        if ap.args.output:
            _generate_config_files(_file, _id, ap.args.threads, ap.args.path)


def _generate_config_files(_file_name: str, _replace_string: str, _threads: int, _outdir: str):
    _config_directory = os.path.join(os.path.dirname(__file__), "config")
    for _config_file in os.listdir(_config_directory):
        _new_file = os.path.join(_outdir, os.path.basename(_config_file))
        cp[os.path.abspath(_config_file), _new_file]()
        _print_and_run(sed["-i", "s/\/path\/to\/%s/\/path\/to\/%s" % (_file_name, _replace_string), _new_file])


def _create_taxonomy_info(mmseqs_db_path: str, outfile: str):
    output_p = open(outfile, "w")
    mmseqs_input_fp = open(mmseqs_db_path + ".lookup", "r")
    for line in mmseqs_input_fp:
        line = line.split()
        output_p.write(line[1] + "\t" + line[1].split("_")[0] + "\n")
    output_p.close()


def _create_tax_db(_out: str, out_dir: str):
    # Create mmseqs input file
    _create_taxonomy_info(_out, os.path.join(out_dir, "mmseqs.input"))
    # Download tax info
    wget["ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", "-O", os.path.join(out_dir, "taxdump.tar.gz")]()
    tar["xzvf", "taxdump.tar.gz"]()
    _print_and_run(
        mmseqs[
            "createtaxdb",
            _out,
            "tmp",
            "--ncbi-tax-dump", out_dir,
            "--tax-mapping-file", os.path.join(out_dir, "mmseqs.input")
        ]
    )


def main(ap: ArgParse, storage_dir: str):
    run(ap, storage_dir)


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("path",),
             {"help": "Download path"}),
            (("-b", "--build"),
             {"help": "Generate required MMseqs2 databases and linear indices, default True", "default": True}),
            (("-x", "--index"),
             {"help": "Generate search index (recommended, but takes a lot of space), default False",
              "default": False}),
            (("-o", "--output"),
             {"help": "Output default config files with included download paths, default True", "default": True}),
            (("-r", "--rewrite"),
             {"help": "Rewrite existing directory, default False", "default": False}),
            (("-t", "--threads"),
             {"help": "Number of threads to use in database generation, default 1", "default": "1"}),
            (("-m", "--max_mem"),
             {"help": "Split memory limit for database generation, default 8G", "default": "8G"}),
        ),
        description="Download required data and build MMseqs2 databases"
    )
    main(_parse_args(_ap), _ap.args.path)
