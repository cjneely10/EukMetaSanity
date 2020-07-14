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


def _parse_args(ap: ArgParse):
    assert os.path.exists(ap.args.path)
    if ap.args.index:
        assert ap.args.build, "Build must be set in order to index!"
    try:
        ap.args.threads = int(ap.args.threads)
    except ValueError as e:
        print(e)
        exit(1)
    return ap


# Download all data with wget from the data.py script - Part of API
def run(ap: ArgParse, pm: PathManager):
    pm.add_dirs(ap.args.path)
    # Download each URL to folder
    for _id, url in data_urls():
        # Download
        pm.add_dirs(ap.args.path, [_id])
        _file = os.path.join(pm.get_dir(ap.args.path), os.path.basename(url.url))
        if not os.path.exists(_file) or ap.args.rewrite:
            wget[url.url, "-O", _file]()
            # Tar/gunzip
            if url.tar:
                tar[url.flags, _file]()
            elif url.gz:
                gunzip[_file]()
            # Generate MMseqs2 database
            _out = os.path.splitext(_file)[0] + "_db"
            if ap.args.build:
                mmseqs[
                    "createdb",
                    _file,
                    _out
                ]()
            # Create MMseqs2 index file
            if ap.args.index:
                mmseqs[
                    "createindex",
                    _out,
                    os.path.join(os.path.basename(_out), "tmp"),
                    "--threads", str(ap.args.threads),
                    "--split-memory-limit", ap.args.max_mem,
                ]()
                mmseqs[
                    "createlinindex",
                    _out,
                    os.path.join(os.path.basename(_out), "tmp"),
                    "--threads", str(ap.args.threads),
                    "--split-memory-limit", ap.args.max_mem,
                ]()
            _generate_config_files(_file, _id, ap.args.threads, ap, pm)


def _generate_config_files(_file_name: str, _replace_string: str, _threads: int, ap: ArgParse, pm: PathManager):
    _config_directory = os.path.join(os.path.dirname(__file__), "config")
    for _config_file in os.listdir(_config_directory):
        _new_file = os.path.join(pm.get_dir(ap.args.path), os.path.basename(_config_file))
        cp[os.path.abspath(_config_file), _new_file]()
        sed["-i", "s/\/path\/to\/%s/\/path\/to\/%s" % (_file_name, _replace_string), _new_file]()


def main(ap: ArgParse, pm: PathManager):
    run(ap, pm)


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
    main(_parse_args(_ap), PathManager(_ap.args.path))
