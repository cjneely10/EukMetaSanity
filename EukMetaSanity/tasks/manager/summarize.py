import os
from plumbum import local
from EukMetaSanity import Task, TaskList


class SummarizeIter(TaskList):
    class Summarize(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = self.input

        def run(self):
            super().run()

    def __init__(self, *args, **kwargs):
        super().__init__(SummarizeIter.Summarize, "summarize", *args, **kwargs)

    def summarize(self, _final_output_dir: str):
        # Create softlinks of final files to output directory
        # Write summary file of reults
        # Store output paths as file for easy loading
        _output = self.output()
        _output_files_list = _output[1]
        _files_prefixes = _output[3]
        if not os.path.exists(_final_output_dir):
            os.makedirs(_final_output_dir)
        _paths_output_file = open(os.path.join(_final_output_dir, "paths_summary.tsv"), "w")
        for _files, _file_prefix in zip(_output_files_list, _files_prefixes):
            _sub_out = os.path.join(_final_output_dir, _file_prefix)
            if not os.path.exists(_sub_out):
                os.makedirs(_sub_out)
            # Copy results to results dir for easier access
            for _file in _files:
                # Write info to file
                if isinstance(_file, dict):
                    sorted_keys = sorted(list(_file.keys()))
                    # Header
                    _paths_output_file.write("".join(("\t".join(["ID"] + sorted_keys), "\n")))
                    # Path info
                    _paths_output_file.write(
                        "".join((
                            "\t".join((
                                _file_prefix,  # Name of record
                                *(os.path.join(_sub_out, str(_file[_f])) for _f in sorted_keys)  # Files for record
                            )), "\n"
                        ))
                    )
                # Generate link of path
                elif isinstance(_file, str):
                    local["ln"]["-srf", _file, _sub_out]()
        _paths_output_file.close()


if __name__ == "__main__":
    pass
