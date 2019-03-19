from abc import ABC, abstractmethod

class AbstractJobScheduler(ABC):
    """
    Abstract class wrapping around a system's job scheduler (e.g. SGE, LSF).
    Defines the minimal set of methods BEERS requires to interact with the
    scheduler to submit, monitor, and kill jobs.
    """

    @staticmethod
    @abstractmethod
    def check_job_status(job_id):
        """
        Determine a job's current run status based on the system's scheduler.

        Parameters
        ----------
        job_id : string
            Identifier used by system's scheduler to uniquely identify the job.

        Returns
        -------
        string
            One of the following:
                RUNNING - job submitted to scheduler and is actively running.
                PENDING - job submitted to scheduler and is pending.
                FAILED - job finished with error status.
                COMPLETED - job finished without error status.

        """
        pass

    @staticmethod
    @abstractmethod
    def submit_job():
        """
        Execute command(s) to submit job to system's scheduler.

        Returns
        -------
        string
            Unique identifier for the submitted job assigned by the system's scheduler.
            Empty string indicates job submission failed.

        """
        pass

    @staticmethod
    @abstractmethod
    def kill_job(job_id):
        """
        Execute command(s) to kill a running job.

        Parameters
        ----------
        job_id : type
            Identifier used by system's scheduler to uniquely identify the job.

        Returns
        -------
        boolean
            True  - Job kill commands executed successfully.
            False - Job kill commands exited with error status.

        """
        pass
