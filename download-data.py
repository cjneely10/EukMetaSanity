import os
from plumbum import local
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.tasks.manager.data import data_urls
from EukMetaSanity.utils.path_manager import PathManager

# Dependencies
mmseqs = local["mmseqs"]
wget = local["wget"]


def _parse_args(ap: ArgParse):
    assert os.path.exists(ap.args.path)
    try:
        ap.args.threads = int(ap.args.threads)
    except ValueError as e:
        print(e)
        exit(1)
    return ap


def _download_fasta_data(ap: ArgParse, pm: PathManager):
    for key, val in data_urls():
        pass


def main(ap: ArgParse, pm: PathManager):
    _download_fasta_data(ap, pm)


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
            (("-m", "--max-mem"),
             {"help": "Split memory limit for database generation, default 8G", "default": "8G"}),
        ),
        description="Download required data and build MMseqs2 databases"
    )
    main(_parse_args(_ap), PathManager(_ap.args.path))
