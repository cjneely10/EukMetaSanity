# Gather all files to parse that match user-passed extensions
def _files_iter(ap: ArgParse, storage_dir: str) -> Generator[str, ArgParse, None]:
    for file in os.listdir(ap.args.fasta_directory):
        for ext in ap.args.extensions:
            if file.endswith(ext):
                yield _simplify_fasta(ap, file, storage_dir)


def _simplify_fasta(ap: ArgParse, file, storage_dir: str) -> str:
    # Simplify FASTA of complex-named sequences
    fasta_file = str(Path(os.path.join(ap.args.fasta_directory, file)).resolve())
    out_file = os.path.join(storage_dir, os.path.basename(os.path.splitext(fasta_file)[0]) + ".fna")
    if os.path.exists(out_file):
        return out_file
    SeqIO.write(SeqIO.parse(fasta_file, "fasta"), out_file, "fasta")
    return out_file


# Get program-needed list of files for this step in pipeline
def _get_list_of_files(summary_file: str) -> Tuple[List[str], List[Dict[str, Dict[str, object]]]]:
    data = json.load(open(summary_file, "r"))
    out_ids = sorted(list(data.keys()))
    out_dict_list = []
    for _id in out_ids:
        to_add = {"root": {}}
        for key, val in data[_id].items():
            if isinstance(val, dict):
                to_add["root"][key] = val
            else:
                to_add[key] = val
        out_dict_list.append(to_add)
    return out_ids, out_dict_list


class InputManager:
    def __init__(self):
        pass
