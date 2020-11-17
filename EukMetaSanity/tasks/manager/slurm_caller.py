import os
from time import sleep
from typing import List, Dict
from plumbum.machines.local import LocalCommand, LocalMachine


class SLURMCaller:

    OUTPUT_SCRIPTS = "slurm-runner.sh"
    FAILED_ID = "failed-job-id"

    def __init__(self, user_id: str, wdir: str, threads: str, cmd: LocalCommand, config_data: Dict[str, str],
                 _local: LocalMachine, slurm_flags: List[str]):
        """

        :param user_id: User id for running slurm job
        :param wdir: Working directory to write slurm script
        :param threads: Number of threads/cpus for job
        :param cmd: Command to run within script
        :param config_data: Dict of data to parse for sbatch script header
        :param _local: local object from which to get sbatch command object
        :param slurm_flags: Added flags to include in script, if any
        """
        self.user_id = user_id
        self.wdir = wdir
        self.threads = threads
        self.cmd = cmd
        self.config = config_data
        self.local = _local
        self.added_flags = slurm_flags

        # Generated job id
        self.job_id: str = SLURMCaller.FAILED_ID
        self.running = False
        # Path of script that is running
        self.script = str(os.path.join(self.wdir, SLURMCaller.OUTPUT_SCRIPTS))
        # Create slurm script in working directory
        self.generate_script()

    def __str__(self):
        return open(self.script, "r").read()

    def __repr__(self):
        return self.__str__()

    # Run script written
    def launch_script(self):
        log_line = str(self.local["sbatch"][self.script]()).split()
        if len(log_line) == 0:
            return
        self.running = self.has_launched(log_line[-1])

    # Call squeue using user id and check if created job id is present
    def is_running(self) -> bool:
        return self.running and self.job_id in str(self.local["squeue"]["-u", self.user_id]).rstrip("\r\n")

    # Parse output from sbatch to see if job id was adequately created
    def has_launched(self, log_line: str) -> bool:
        try:
            self.job_id = str(int(log_line))
        except ValueError:
            return False
        return True

    # check if a particular job is still running
    def has_ended(self) -> bool:
        return not self.is_running()

    # Create slurm script to run program
    def generate_script(self):
        # Initialize file
        fp = open(self.script, "w")
        fp.write("#!/bin/bash\n\n")

        # Write all header lines
        fp.write(SLURMCaller.create_header_line("--nodes", "1"))
        fp.write(SLURMCaller.create_header_line("--tasks", "1"))
        fp.write(SLURMCaller.create_header_line("--cpus-per-task", self.threads))
        fp.write(SLURMCaller.create_header_line("--mem", self.config["MEMORY"]))
        fp.write(SLURMCaller.create_header_line("--time", self.config["TIME"]))
        # Write additional header lines passed in by user
        for added_arg in self.added_flags:
            fp.write(SLURMCaller.create_header_line(*added_arg.split("=")))
        fp.write("\n")
        # Write command to run
        fp.write("".join((str(self.cmd), "\n")))
        fp.close()

    @staticmethod
    def create_header_line(param: str, arg: str) -> str:
        return "#SBATCH %s=%s\n" % (param, arg)

    # Call will run script using sbatch and periodically check for when to stop running
    def __call__(self, *args, **kwargs):
        self.launch_script()
        while self.is_running():
            sleep(600)  # Wait 10 minutes in between checking if still running
