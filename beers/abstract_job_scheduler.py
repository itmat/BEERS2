from abc import ABC, abstractmethod

class AbstractJobScheduler(ABC):
    """
    Abstract class wrapping around a system's job scheduler (e.g. SGE, LSF).
    Defines the minimal set of methods BEERS requires to interact with the
    scheduler to submit, monitor, and kill jobs.
    """

    @staticmethod
    @abstractmethod
    def check_job_status(job_id, additional_args):
        """
        Execute command(s) to determine a job's current run status based on the
        system's scheduler.

        Parameters
        ----------
        job_id : string
            Identifier used by system's scheduler to uniquely identify the job.
        additional_args : string
            Additonal arguments to provide to the command(s) used to check the
            job's status.

        Returns
        -------
        string
            One of the following:
                RUNNING - according to scheduler and the job is actively running.
                PENDING - according to scheduler and the job is pending.
                FAILED - according to scheduler the job finished with error status.
                COMPLETED - according to scheduler the job finished without error status.
                ERROR - could not retrieve job status from scheduler.

        """
        pass

    @staticmethod
    @abstractmethod
    def submit_job(job_command, job_name, stdout_logfile, stderr_logfile,
                   memory_in_mb, num_processors, additional_args):
        """
        Execute command(s) to submit job to system's scheduler.

        Parameters
        ----------
        job_command : string
            Full command to execute job when run from the command line.
        job_name : string
            Name assigned to job by scheduler.
        stdout_logfile : string
            Full path to file where stdout from the scheduler/command should be
            stored.
        stderr_logfile : string
            Full path to file where stderr from the scheduler/command should be
            stored.
        memory_in_mb : int
            Memory (in Mb) to request for running the job.
        num_processors : int
            Number of processors/cores to request for running the job.
        additional_args : string
            Additonal arguments to provide to the command(s) used to submit the job.

        Returns
        -------
        string
            Unique identifier for the submitted job assigned by the system's scheduler.
            "ERROR" string indicates job submission failed.

        """
        pass

    @staticmethod
    @abstractmethod
    def kill_job(job_id, additional_args):
        """
        Execute command(s) to kill a running job.

        Parameters
        ----------
        job_id : type
            Identifier used by system's scheduler to uniquely identify the job.
        additional_args : string
            Additonal arguments to provide to the command(s) used to kill the job.

        Returns
        -------
        boolean
            True  - Job kill commands executed successfully.
            False - Job kill commands exited with error status.

        """
        pass
