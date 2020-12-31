"""
Module holds logic for running a dask distributed task within a SLURM job
"""

import os
from time import sleep
from typing import List, Tuple
from plumbum.machines.local import LocalCommand, LocalMachine


class SLURMCaller:
    """ SLURMCaller handles running a program on a SLURM cluster

    Overrides __str__ and __call__ method to allow script command to be displayed within config
    and called by dask task manager, respectively.

    """
    OUTPUT_SCRIPTS = "slurm-runner.sh"
    FAILED_ID = "failed-job-id"

    def __init__(self, user_id: str, wdir: str, threads: str, cmd: LocalCommand,
                 memory: str, time: str, _local: LocalMachine, slurm_flags: List[Tuple[str, str]]):
        """ Generate SLURMCaller object using user metadata gathered from SLURM config section and
        the task's own metadata

        :param user_id: user-id argument in SLURM template
        :param wdir: Job working directory
        :param threads: Total number of cpus to give to job
        :param cmd: Command to run in SLURM template
        :param memory: Amount of memory to request for job
        :param time: Max allowed time to run job
        :param _local: Reference to LocalMachine object to use to run command
        :param slurm_flags: SLURM-specific flags to pass to top of SLURM job script (e.g. SBATCH flags)
        """
        self.user_id = user_id
        self.wdir = wdir
        self.threads = threads
        self.cmd = cmd
        self.memory = memory
        self.time = time
        self.local = _local
        self.added_flags = slurm_flags

        # Generated job id
        self.job_id: str = SLURMCaller.FAILED_ID
        # Initialize as not running
        self.running = False
        # Path of script that will run
        self.script = str(os.path.join(self.wdir, SLURMCaller.OUTPUT_SCRIPTS))
        # Create slurm script in working directory
        self._generate_script()

    def __str__(self) -> str:
        """ Return contents of script as a string

        :return: Contents of SLURM script
        """
        return open(self.script, "r").read()

    def __repr__(self) -> str:
        """ REPL version of string

        :return: Contents of SLURM script
        """
        return self.__str__()

    def _launch_script(self):
        """ Run generated script using sbatch

        Check if loaded properly and store in object data if properly launched
        """
        # Call script using sbatch
        log_line = str(self.local["sbatch"][self.script]()).split()
        # If no output to stdout/err, return
        if len(log_line) == 0:
            return
        # Otherwise set running status based on contents of stdout/err
        self.running = self._has_launched(log_line[-1])

    # Call squeue using user id and check if created job id is present
    def _is_running(self) -> bool:
        """ Method will check if the job_id created by launching the sbatch script is still
        visible within a user's `squeue`.

        :return: Task still running (true) or has completed/failed to start (false)
        """
        return self.running and self.job_id in str(self.local["squeue"]["-u", self.user_id]())

    def _has_launched(self, log_line: str) -> bool:
        """ Parse output from sbatch to see if job id was adequately created

        :param log_line: Output from running sbatch that should be a parsable integer
        :return: True if launched properly, or False if unable to parse/task was not launched properly
        """
        try:
            self.job_id = str(int(log_line))
        except ValueError:
            return False
        return True

    def _generate_script(self):
        """ Create script using passed command within SLURM format

        Include all additional flags to SLURM

        Flags to command should already be passed to command itself
        """
        # Initialize file
        file_ptr = open(self.script, "w")
        file_ptr.write("#!/bin/bash\n\n")

        # Write all header lines
        file_ptr.write(SLURMCaller._create_header_line("--nodes", "1"))
        file_ptr.write(SLURMCaller._create_header_line("--tasks", "1"))
        file_ptr.write(SLURMCaller._create_header_line("--cpus-per-task", self.threads))
        file_ptr.write(SLURMCaller._create_header_line("--mem", self.memory))
        file_ptr.write(SLURMCaller._create_header_line("--time", self.time))
        # Write additional header lines passed in by user
        for added_arg in self.added_flags:
            file_ptr.write(SLURMCaller._create_header_line(*added_arg))
        file_ptr.write("\n")
        # Write command to run
        file_ptr.write("".join((str(self.cmd), "\n")))
        file_ptr.close()

    @staticmethod
    def _create_header_line(param: str, arg: str) -> str:
        """ Creates a header line in a SLURM script of form SBATCH key=val

        :param param: key
        :param arg: val
        :return: str of header line (with newline added)
        """
        return "#SBATCH %s=%s\n" % (param, arg)

    def __call__(self, *args, **kwargs):
        """ Call will run script using sbatch and periodically check for when to stop running

        :param args: Any args passed
        :param kwargs: Any kwargs passed
        """
        self._launch_script()
        while self._is_running():
            sleep(60)  # Wait 1 minute in between checking if still running
