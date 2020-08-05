import os
from EukMetaSanity import Task, TaskList, program_catch

"""
Call MAKER3 pipeline

"""


class MakerIter(TaskList):
    class Maker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".all.maker.gff"),  # MAKER output gff
                os.path.join(self.wdir, self.record_id + ".proteins.fasta"),  # Proteins output
                os.path.join(self.wdir, self.record_id + ".transcripts.fasta"),  # CDS output
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.reformat_repeats_gff3()
            # # Run initial MAKER
            self.generate_ctl_file()
            # Run maker using edited config files
            self.log_and_run(
                self.program_mpi[
                    "-n", self.threads, self.program,
                    "-base", os.path.join(self.wdir, self.record_id),
                    os.path.join(self.wdir, "maker_opts.ctl"),
                    os.path.join(self.wdir, "maker_bopts.ctl"),
                    os.path.join(self.wdir, "maker_exe.ctl"),
                ] | self.local["tee"][os.path.join(self.wdir, "maker.log")]
            )
            self.merge_output()
            # # Run on remaining transcriptomes that were assembled in assemble task

        # Generate either using 'maker' or 'initial' type
        def generate_ctl_file(self, file_path: str, file_type: str = "initial"):
            # Create base file and move to wdir
            if not os.path.exists("maker_opts.ctl"):
                self.log_and_run(self.program["-CTL"])
                self.log_and_run(self.local["cp"]["*.ctl", self.wdir])
            opts_file = os.path.join(self.wdir, "maker_opts.ctl")
            # Add proper paths to file
            # Edit path to genome
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("genome=", "genome=%s" % self.input[4]),
                    opts_file,
                ]
            )
            # Incorporate search for simple repeats with MAKER, use existing models from Run
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("model_org=", "model_org=simple"),
                    opts_file,
                ]
            )
            self.log_and_run(
                self.local["sed"][
                    "-i",
                    "s/%s/%s/" % (
                        "rm_gff=",
                        "rm_gff=%s" % os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
                    ),
                    opts_file,
                ]
            )
            # Add nr evidence
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("other_gff=", "other_gff=%s" % self.input[2]),
                    opts_file,
                ]
            )
            # Parse additional evidence as needed
            prot_file = self.get_data_file(self.proteins)
            if prot_file is not None:
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/%s/%s/" % ("protein=", "protein=%s" % prot_file),
                        opts_file
                    ]
                )
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/protein2genome=0/protein2genome=1/",
                        opts_file
                    ]
                )
            # Integrate assembed transcriptomes, if available
            est_file = self.get_data_file(self.est)
            if est_file is not None:
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/%s/%s/g" % ("est=", "est=%s" % est_file),
                        opts_file
                    ]
                )
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/est2genome=0/est2genome=1/",
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

        def merge_output(self):
            self.log_and_run(
                self.program_gff3merge[
                    "-s", "-d", os.path.join(self.wdir, self.record_id + "_master_datastore_index.log")
                ] > os.path.join(self.wdir, self.record_id + ".all.maker.gff")
            )
            self.log_and_run(
                self.program_fastamerge[
                    "-d", os.path.join(self.wdir, self.record_id + "_master_datastore_index.log")
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(MakerIter.Maker, "maker", *args, **kwargs)


if __name__ == "__main__":
    pass
