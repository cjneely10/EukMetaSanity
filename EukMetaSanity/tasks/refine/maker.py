import os
from typing import Tuple, List
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import prefix

"""
Call MAKER3 pipeline

"""


class MakerIter(TaskList):
    class Maker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            out = []
            for i, assembled_transcriptome in enumerate(self.input[-1]):
                run_id = os.path.join(self.wdir, prefix(assembled_transcriptome) + self.record_id)
                for ext in (".all.maker.gff", ".proteins.fasta", ".transcripts.fasta"):
                    out.append(os.path.join(self.wdir, run_id + ext))
            self.output = [
                *self.input,
                *out,  # MAKER output paths
                self.input[-2],  # Forward the list of reads
                out,  # List of MAKER results files
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.reformat_repeats_gff3()
            # Run MAKER on assembled transcriptomes
            out_files = None
            prot_file = self.get_data_file(self.proteins)
            if prot_file is None:
                prot_file = ""
            for i, assembled_transcriptome in enumerate(self.input[-1]):
                out_files = self.maker(
                    *self.generate_ctl_file(
                        [
                            ("genome", self.input[4]),
                            ("model_org", ("simple" if i == 0 else "")),
                            ("rm_gff", (
                                os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
                                if i == 0 else out_files[2]
                            )),
                            ("other_gff", (self.input[2] if i == 0 else out_files[0])),  # Use est evidence from prior
                            *(("protein", prot_file) if i == 0 else out_files[1]),
                            ("est", assembled_transcriptome),
                        ],
                        i + 1
                    ),
                    os.path.join(self.wdir, prefix(assembled_transcriptome) + self.record_id)
                )
            # # Run on remaining transcriptomes that were assembled in assemble task

        @program_catch
        def maker(self, maker_opts: str, maker_bopts: str, maker_exe: str, run_id: str) -> List[str]:
            # Run maker using edited config files
            self.log_and_run(
                self.program_mpi[
                    "-n", self.threads, self.program,
                    "-base", run_id,
                    maker_opts,
                    maker_bopts,
                    maker_exe,
                ] | self.local["tee"][os.path.join(self.wdir, "maker-%s.log" % run_id)]
            )
            return self.merge_output(run_id)

        # Generate either using 'maker' or 'initial' type
        def generate_ctl_file(self, replace_tuples: List[Tuple[str, str]], _round: int) -> Tuple[str, str, str]:
            # Create base file and move to wdir
            if not os.path.exists("maker_opts.ctl"):
                self.log_and_run(self.program["-CTL"])
            self.log_and_run(self.local["cp"]["*.ctl", self.wdir])
            # Rename files based on run number
            os.replace(os.path.join(self.wdir, "maker_opts.ctl"), os.path.join(self.wdir, "maker_opts%i.ctl" % _round))
            os.replace(os.path.join(self.wdir, "maker_exe.ctl"), os.path.join(self.wdir, "maker_exe%i.ctl" % _round))
            os.replace(
                os.path.join(self.wdir, "maker_bopts.ctl"), os.path.join(self.wdir, "maker_bopts%i.ctl" % _round)
            )
            opts_file = os.path.join(self.wdir, "maker_opts%i.ctl" % _round)
            # Replace value in config file with valid concatenated key=path value
            for replace_tuple in replace_tuples:
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/%s/%s/" % ("%s=" % replace_tuple[0], "%s=%s" % (replace_tuple[0], replace_tuple[1])),
                        opts_file,
                    ]
                )
                if replace_tuple[0] in ("protein", "est"):
                    self.log_and_run(
                        self.local["sed"][
                            "-i", "s/%s2genome=0/%s2genome=1/" % (replace_tuple[0], replace_tuple[0]),
                            opts_file
                        ]
                    )
            # Parse user args into config file
            assert len(self.added_flags) // 2 == 0
            for _i in range(0, len(self.added_flags) - 1, 2):
                self.log_and_run(
                    self.local["sed"][
                        "-i",
                        "s/%s/%s" % (
                            self.added_flags[_i] + "=",
                            self.added_flags[_i] + "=" + self.added_flags[_i + 1]
                        ),
                        opts_file
                    ]
                )
            return (os.path.join(self.wdir, "maker_opts%i.ctl" % _round),
                    os.path.join(self.wdir, "maker_bopts%i.ctl" % _round),
                    os.path.join(self.wdir, "maker_exe%i.ctl" % _round))

        def reformat_repeats_gff3(self):
            # Isolate complex repeats
            self.log_and_run(
                self.local["grep"][
                    "-v", "-e", "Satellite", "-e", ")n", "-e", "-rich",
                    self.input[3]
                ] > os.path.join(self.wdir, self.record_id + ".mask.complex.gff3")
            )
            # Reformat to work with MAKER
            self.log_and_run(
                self.local["cat"][os.path.join(self.wdir, self.record_id + ".mask.complex.gff3")] |
                self.local["perl"][
                    "-ane",
                    r'$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", '
                    r'@F)."\n"} print $_',
                ] > os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
            )

        def get_data_file(self, file_path: str):
            if not os.path.exists(file_path):
                return
            with open(file_path, "r") as w:
                for line in w:
                    if line.startswith(self.record_id):
                        return line.rstrip("\r\n").split("\t")[1]

        def merge_output(self, run_id: str) -> List[str]:
            # Generate gff3 file
            self.log_and_run(
                self.program_gff3merge[
                    "-n", "-s", "-d", os.path.join(self.wdir, run_id + "_master_datastore_index.log")
                ] > os.path.join(self.wdir, run_id + ".all.maker.gff")
            )
            # Separate based on evidence type
            out_files = []
            for val in ("est2genome", "protein2genome", "repeat"):
                _file = os.path.join(self.wdir, run_id + ".%s.maker.gff" % val)
                self.log_and_run(
                    self.local["awk"][
                        '{ if ($2 == "%s") print $0 }' % val,
                        os.path.join(self.wdir, run_id + ".all.maker.gff")
                    ] > _file
                )
                if os.path.exists(_file):
                    if os.path.getsize(_file) == 0:
                        self.local["rm"][_file]()
                        out_files.append("")
                    else:
                        out_files.append(_file)
            # Generate FASTA files
            self.log_and_run(
                self.program_fastamerge[
                    "-d", os.path.join(self.wdir, run_id + "_master_datastore_index.log")
                ]
            )
            return out_files

    def __init__(self, *args, **kwargs):
        super().__init__(MakerIter.Maker, "maker", *args, **kwargs)


if __name__ == "__main__":
    pass
