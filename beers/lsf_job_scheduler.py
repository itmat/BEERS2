import re
import subprocess
from beers.abstract_job_scheduler import AbstractJobScheduler

class LsfJobScheduler(AbstractJobScheduler):
    """
    Wrapper around the Load Sharing Facility (LSF) scheduler. Provides methods
    for monitoring job status, as well as submitting and killing jobs.
    """

    #Regular expression for parsing bjobs output (including header)
    LSF_BJOBS_OUTPUT_PATTERN = re.compile(r'''JOBID\s+USER\s+STAT\s+QUEUE\s+FROM_HOST\s+EXEC_HOST\s+JOB_NAME\s+SUBMIT_TIME\n(?P<job_id>\d+?)\s+\S+\s+(?P<job_status>\S+?)\s+.*''')

    #Default command used to submit jobs.
    DEFAULT_BSUB_COMMAND = ('bsub -J {job_name}'
                            ' -n {num_processors}'
                            ' -R \"span[hosts=1]\"'
                            ' -M {mem_usage_in_mb}'
                            ' -R \"rusage[mem={mem_usage_in_mb}]\"'
                            ' -oo {stdout_logfile}'
                            ' -eo {stderr_logfile}')

    @staticmethod
    def check_job_status(job_id):
        """
        Return status of given job in the LSF queue. This operation performed
        using the "bjobs" command.

        Parameters
        ----------
        job_id : string
            Unique LSF job id.

        Returns
        -------
        string
            One of the following:
                RUNNING - according to LSF scheduler and the job is actively running.
                PENDING - according to LSF scheduler and the job is pending.
                FAILED - according to LSF scheduler the job finished with error status.
                COMPLETED - according to LSF scheduler the job finished without error status.
                ERROR - could not retrieve job status from LSF scheduler.

        """

        job_status = "ERROR"

        bjobs_result = subprocess.run(' '.join([f"bjobs {job_id}"]), shell=True, check=True, stdout=subprocess.PIPE, encoding="ascii")

        #Here's some code just using string split to try to get job status
        #Skip first line bjobs output, since it just contains the header info.
        #job_status = bjobs_result.stdout.split("\n")[1].split()[2]

        if LsfJobScheduler.LSF_BJOBS_OUTPUT_PATTERN.match(bjobs_result.stdout):

            lsf_job_status = LsfJobScheduler.LSF_BJOBS_OUTPUT_PATTERN.match(bjobs_result.stdout).group("job_status")

            if lsf_job_status == "RUN":
                job_status = "RUNNING"
            elif lsf_job_status == "PEND" or lsf_job_status == "WAIT":
                job_status = "PENDING"
            elif lsf_job_status == "EXIT":
                job_status = "FAILED"
            elif lsf_job_status == "DONE":
                job_status = "COMPLETED"

        return job_status

    @staticmethod
    def submit_job(job_command, submission_args=None):
        """
        Submit given job using the bsub command and return ID assigned to job by
        the LSF scheduler.

        Parameters
        ----------
        job_command : string
            Full command to execute job when run from the command line. It cannot
            contain any unix output redirection (i.e. useing ">" or "2>") unless
            the entire command is enclosed in single-quotes.
        submission_args : dict
            Arguments and corresponding values to pass to the bsub command. These
            will override any default arguments used in the bsub command.

        Returns
        -------
        string
            Unique identifier for the submitted job assigned by the LSF scheduler.
            Empty string indicates job submission failed.

        """
        pass

    @staticmethod
    def kill_job(job_id, kill_args=None):
        """
        Kill given job using the bkill command.

        Parameters
        ----------
        job_id : string
            Unique LSF job id.
        kill_args : string
            Additional arguments to pass to bkill command.

        Returns
        -------
        boolean
            True  - bkill command executed successfully.
            False - bkill command exited with error status.

        """
        pass
