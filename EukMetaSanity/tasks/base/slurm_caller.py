import os
from time import sleep
from typing import List, Dict, Tuple, Optional, Union
from plumbum.machines.local import LocalCommand, LocalMachine
from EukMetaSanity.tasks.base.config_manager import ConfigManager


class SLURMCaller:
    """ SLURMCaller handles running a program on a SLURM cluster

    Overrides __str__ and __call__ method to allow script command to be displayed within config
    and called by dask task manager, respectively.

    """
    OUTPUT_SCRIPTS = "slurm-runner.sh"
    FAILED_ID = "failed-job-id"

    def __init__(self, user_id: str, wdir: str, threads: str, cmd: LocalCommand,
                 config_data: Dict[str, Dict[str, Union[str, dict]]], _local: LocalMachine,
                 slurm_flags: List[Tuple[str, str]], time_override: Optional[str]):
        self.user_id = user_id
        self.wdir = wdir
        self.threads = threads
        self.cmd = cmd
        self.config = config_data
        self.local = _local
        self.added_flags = slurm_flags
        self.time_override = time_override

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
        """ Run generated script using sbatch

        Check if loaded properly and store in object data if properly launched
        """
        log_line = str(self.local["sbatch"][self.script]()).split()
        if len(log_line) == 0:
            return
        self.running = self.has_launched(log_line[-1])

    # Call squeue using user id and check if created job id is present
    def is_running(self) -> bool:
        """ Method will check if the job_id created by launching the sbatch script is still
        visible within a user's `squeue`.

        :return: Task still running (true) or has completed/failed to start (false)
        """
        return self.running and self.job_id in str(self.local["squeue"]["-u", self.user_id]())

    def has_launched(self, log_line: str) -> bool:
        """ Parse output from sbatch to see if job id was adequately created

        :param log_line: Output from running sbatch that should be a parsable integer
        :return: True if launched properly, or False if unable to parse/task was not launched properly
        """
        try:
            self.job_id = str(int(log_line))
        except ValueError:
            return False
        return True

    # check if a particular job is still running
    def has_ended(self) -> bool:
        """ Checks if a task is still running on a user's `squeue`

        :return: True is job stopped running, or False if still running
        """
        return not self.is_running()

    # Create slurm script to run program
    def generate_script(self):
        """ Create script using passed command within SLURM format

        Include all additional flags to SLURM

        Flags to command should already be passed to command itself
        """
        # Initialize file
        fp = open(self.script, "w")
        fp.write("#!/bin/bash\n\n")

        # Write all header lines
        fp.write(SLURMCaller.create_header_line("--nodes", "1"))
        fp.write(SLURMCaller.create_header_line("--tasks", "1"))
        fp.write(SLURMCaller.create_header_line("--cpus-per-task", self.threads))
        fp.write(SLURMCaller.create_header_line("--mem", str(self.config[ConfigManager.MEMORY])))
        if ConfigManager.TIME in self.config.keys():
            if self.time_override is not None:
                fp.write(SLURMCaller.create_header_line("--time", self.time_override))
            else:
                fp.write(SLURMCaller.create_header_line("--time", str(self.config[ConfigManager.TIME])))
        # Write additional header lines passed in by user
        for added_arg in self.added_flags:
            fp.write(SLURMCaller.create_header_line(*added_arg))
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
            sleep(60)  # Wait 1 minute in between checking if still running
