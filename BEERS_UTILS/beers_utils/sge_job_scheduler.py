import re
import subprocess
from beers_utils.abstract_job_scheduler import AbstractJobScheduler

class SgeJobScheduler(AbstractJobScheduler):
    """
    Wrapper around the Sun Grid Engine (SGE) scheduler. Provides methods
    for monitoring job status, as well as submitting and killing jobs.
    """

    #Regular expression for parsing qstat output (including header)
    SGE_QSTAT_OUTPUT_PATTERN = re.compile(r'''job-ID\s+prior\s+name\s+user\s+state\s+submit/start at\s+queue\s+slots\s+ja-task-ID\n(?P<job_id>\d+?)\s+\S+\s+(?P<job_status>\S+?)\s+.*''')

    #Regular expression for recognizing and extracting SGE job IDs following qsub.
    SGE_QSUB_OUTPUT_PATTERN = re.compile(r'Your job (?P<job_id>\d+?) \(.*\)  has been submitted')

    #Regular expression for parsing qdel output.
    SGE_QDEL_OUTPUT_PATTERN = re.compile(r'.* has registered the job (?P<job_id>\d+?) for deletion(?P<qdel_message>[^\n]+?)\n?$')

    #Default command used to check job status.
    DEFAULT_QSTAT_COMMAND = ('qstat {qstat_args} {job_id}')

    #Default command used to submit job.
    DEFAULT_QSUB_COMMAND = ('qsub -N \"{job_name}\"'
                            ' -V -cwd'
                            ' -n {num_processors}'
                            ' -R \"span[hosts=1]\"'
                            ' -l h_vmem={mem_usage_in_mb}M')

    #Default command used to kill job.
    DEFAULT_QDEL_COMMAND = ('qdel {qdel_args} {job_id}')

    @staticmethod
    def check_job_status(job_id, additional_args=""):
        """
        Return status of given job in the SGE queue. This operation performed
        using the "qstat" command.

        Parameters
        ----------
        job_id : string
            Unique SGE job id.
        additional_args : string
            Additional arguments to provide to the qstat command. Default: empty string.

        Returns
        -------
        string
            One of the following:
                RUNNING - according to SGE scheduler and the job is actively running.
                PENDING - according to SGE scheduler and the job is pending.
                FAILED - according to SGE scheduler the job finished with error status.
                COMPLETED - according to SGE scheduler the job finished without error status.
                ERROR - could not retrieve job status from SGE scheduler.

        """

        job_status = "ERROR"

        qstat_command = SgeJobScheduler.DEFAULT_QSTAT_COMMAND.format(job_id=job_id, qstat_args=additional_args)
        qstat_result = subprocess.run(qstat_command, shell=True, check=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding="ascii")

        #Here's some code just using string split to try to get job status
        #Skip first line qstat output, since it just contains the header info.
        #job_status = qstat_result.stdout.split("\n")[1].split()[2]

        if SgeJobScheduler.SGE_QSTAT_OUTPUT_PATTERN.match(qstat_result.stdout):

            sge_job_status = SgeJobScheduler.SGE_QSTAT_OUTPUT_PATTERN.match(qstat_result.stdout).group("job_status")

            if sge_job_status == "r":
                job_status = "RUNNING"
            elif sge_job_status == "qw" or sge_job_status == "WAIT":
                job_status = "PENDING"
            elif sge_job_status == "EXIT" or sge_job_status == "ZOMBI" or sge_job_status == "UNKWN":
                job_status = "FAILED"
            elif sge_job_status == "DONE":
                job_status = "COMPLETED"

        return job_status

    @staticmethod
    def submit_job(job_command, job_name, stdout_logfile=None, stderr_logfile=None,
                   memory_in_mb=6000, num_processors=1, additional_args=""):
        """
        Submit given job using the qsub command and return ID assigned to job by
        the SGE scheduler.

        Parameters
        ----------
        job_command : string
            Full command to execute job when run from the command line. It cannot
            contain any unix output redirection (i.e. useing ">" or "2>") unless
            the entire command is enclosed in single-quotes.
        job_name : string
            Name assigned to job by SGE. Specified with '-J' argument.
        stdout_logfile : string
            Full path to file where SGE/job stdout should be stored. Specified
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
            Arguments and corresponding values to pass to the qsub command.
            Default: empty string.

        Returns
        -------
        string
            Unique identifier for the submitted job assigned by the SGE scheduler.
            "ERROR" string indicates job submission failed.

        """
        job_id = "ERROR"

        qsub_command = SgeJobScheduler.DEFAULT_QSUB_COMMAND.format(job_name=job_name,
                                                                   num_processors=num_processors,
                                                                   mem_usage_in_mb=memory_in_mb)
        if stdout_logfile:
            qsub_command += f" -oo {stdout_logfile}"
        if stderr_logfile:
            qsub_command += f" -eo {stderr_logfile}"

        qsub_result = subprocess.run(' '.join([qsub_command, additional_args, job_command]),
                                     shell=True, check=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, encoding="ascii")
        if SgeJobScheduler.SGE_QSUB_OUTPUT_PATTERN.match(qsub_result.stdout):
            job_id = SgeJobScheduler.SGE_QSUB_OUTPUT_PATTERN.match(qsub_result.stdout).group('job_id')

        return job_id

    @staticmethod
    def kill_job(job_id, additional_args=""):
        """
        Kill given job using the qdel command.

        Parameters
        ----------
        job_id : string
            Unique SGE job id.
        additional_args : string
            Additional arguments to pass to qdel command. Default: empty string.

        Returns
        -------
        boolean
            True  - qdel command executed successfully.
            False - qdel command exited with error status.

        """
        qdel_status = False

        qdel_command = SgeJobScheduler.DEFAULT_QDEL_COMMAND.format(job_id=job_id, qdel_args=additional_args)
        #Note, set check=False here, since qdel command exits with an error code
        #if the given job does not exist or has already completed. These cases
        #are handled by the code afterward, so we don't need to check here (which
        #causes the whole run() statement to throw an exception).
        qdel_result = subprocess.run(qdel_command, shell=True, check=False,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding="ascii")
        qdel_message = SgeJobScheduler.SGE_QDEL_OUTPUT_PATTERN.match(qdel_result.stdout).group('qdel_message')

        if qdel_message == " is being terminated" or qdel_message == ": Job has already finished":
            qdel_status = True

        #TODO: Wait and check job status to make sure job was really killed. Could
        #      escalate to different qdel arguments if the job fails to die, or
        #      simply throw an exception and let the calling code decided how
        #      to handle this.

        return qdel_status
