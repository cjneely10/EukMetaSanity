import os
from plumbum import local
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.utils.path_manager import PathManager

# URLS to gather
ORTHODB_URL = "https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz"

# Dependencies
mmseqs = local["mmseqs"]
wget = local["wget"]


def _parse_args(ap: ArgParse):
    assert os.path.exists(ap.args.path)


def main(ap: ArgParse):
    pass


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
             {"help": "Rewrite existing directory, default False", "default": False})
        ),
        description="Download required data and build MMseqs2 databases"
    )
    _parse_args(_ap)
    main(_ap)
