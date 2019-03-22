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

    #Regular expression for recognizing and extracting lsf job IDs following bsub.
    LSF_BSUB_OUTPUT_PATTERN = re.compile(r'Job <(?P<job_id>\d+?)> is submitted .*')

    #Regular expression for parsing bkill output.
    LSF_BKILL_OUTPUT_PATTERN = re.compile(r'Job <(?P<job_id>\d+?)>(?P<bkill_message>[^\n]+?)\n?$')

    #Default command used to check job status.
    DEFAULT_BJOBS_COMMAND = ('bjobs {bjobs_args} {job_id}')

    #Default command used to submit job.
    DEFAULT_BSUB_COMMAND = ('bsub -J \"{job_name}\"'
                            ' -n {num_processors}'
                            ' -R \"span[hosts=1]\"'
                            ' -M {mem_usage_in_mb}'
                            ' -R \"rusage[mem={mem_usage_in_mb}]\"')

    #Default command used to kill job.
    DEFAULT_BKILL_COMMAND = ('bkill {bkill_args} {job_id}')

    @staticmethod
    def check_job_status(job_id, additional_args=""):
        """
        Return status of given job in the LSF queue. This operation performed
        using the "bjobs" command.

        Parameters
        ----------
        job_id : string
            Unique LSF job id.
        additional_args : string
            Additional arguments to provide to the bjobs command. Default: empty string.

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

        bjobs_command = LsfJobScheduler.DEFAULT_BJOBS_COMMAND.format(job_id=job_id, bjobs_args=additional_args)
        bjobs_result = subprocess.run(bjobs_command, shell=True, check=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding="ascii")

        #Here's some code just using string split to try to get job status
        #Skip first line bjobs output, since it just contains the header info.
        #job_status = bjobs_result.stdout.split("\n")[1].split()[2]

        if LsfJobScheduler.LSF_BJOBS_OUTPUT_PATTERN.match(bjobs_result.stdout):

            lsf_job_status = LsfJobScheduler.LSF_BJOBS_OUTPUT_PATTERN.match(bjobs_result.stdout).group("job_status")

            if lsf_job_status == "RUN":
                job_status = "RUNNING"
            elif lsf_job_status == "PEND" or lsf_job_status == "WAIT":
                job_status = "PENDING"
            elif lsf_job_status == "EXIT" or lsf_job_status == "UNKWN":
                job_status = "FAILED"
            elif lsf_job_status == "DONE":
                job_status = "COMPLETED"

        return job_status

    @staticmethod
    def submit_job(job_command, job_name, stdout_logfile=None, stderr_logfile=None,
                   memory_in_mb=6000, num_processors=1, additional_args=""):
        """
        Submit given job using the bsub command and return ID assigned to job by
        the LSF scheduler.

        Parameters
        ----------
        job_command : string
            Full command to execute job when run from the command line. It cannot
            contain any unix output redirection (i.e. useing ">" or "2>") unless
            the entire command is enclosed in single-quotes.
        job_name : string
            Name assigned to job by LSF. Specified with '-J' argument.
        stdout_logfile : string
            Full path to file where LSF/job stdout should be stored. Specified
            using the "-oo" argument. Default: None.
        stderr_logfile : string
            Full path to file where job stderr should be stored. Specified using
            the "-eo" argument. Default: None.
        memory_in_mb : int
            Memory (in Mb) to request for running the job. Specified using both
            the '-M' and '-R "rusage[mem=]"' arguments. Default: 6000.
        num_processors : int
            Number of processor units to request for running the job. Specified
            using the '-n' argument. Default: 1.
        additional_args : string
            Arguments and corresponding values to pass to the bsub command.
            Default: empty string.

        Returns
        -------
        string
            Unique identifier for the submitted job assigned by the LSF scheduler.
            "ERROR" string indicates job submission failed.

        """
        job_id = "ERROR"

        bsub_command = LsfJobScheduler.DEFAULT_BSUB_COMMAND.format(job_name=job_name,
                                                                   num_processors=num_processors,
                                                                   mem_usage_in_mb=memory_in_mb)
        if stdout_logfile:
            bsub_command += f" -oo {stdout_logfile}"
        if stderr_logfile:
            bsub_command += f" -eo {stderr_logfile}"

        bsub_result = subprocess.run(' '.join([bsub_command, additional_args, job_command]),
                                     shell=True, check=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, encoding="ascii")
        if LsfJobScheduler.LSF_BSUB_OUTPUT_PATTERN.match(bsub_result.stdout):
            job_id = LsfJobScheduler.LSF_BSUB_OUTPUT_PATTERN.match(bsub_result.stdout).group('job_id')

        return job_id

    @staticmethod
    def kill_job(job_id, additional_args=""):
        """
        Kill given job using the bkill command.

        Parameters
        ----------
        job_id : string
            Unique LSF job id.
        additional_args : string
            Additional arguments to pass to bkill command. Default: empty string.

        Returns
        -------
        boolean
            True  - bkill command executed successfully.
            False - bkill command exited with error status.

        """
        bkill_status = False

        bkill_command = LsfJobScheduler.DEFAULT_BKILL_COMMAND.format(job_id=job_id, bkill_args=additional_args)
        #Note, set check=False here, since bkill command exits with an error code
        #if the given job does not exist or has already completed. These cases
        #are handled by the code afterward, so we don't need to check here (which
        #causes the whole run() statement to throw an exception).
        bkill_result = subprocess.run(bkill_command, shell=True, check=False,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding="ascii")
        bkill_message = LsfJobScheduler.LSF_BKILL_OUTPUT_PATTERN.match(bkill_result.stdout).group('bkill_message')

        if bkill_message == " is being terminated" or bkill_message == ": Job has already finished":
            bkill_status = True

        #TODO: Wait and check job status to make sure job was really killed. Could
        #      escalate to different bkill arguments if the job fails to die, or
        #      simply throw an exception and let the calling code decided how
        #      to handle this.

        return bkill_status
