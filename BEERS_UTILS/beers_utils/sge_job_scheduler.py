import re
import time
import subprocess
from beers_utils.abstract_job_scheduler import AbstractJobScheduler

# TODO: Consider shifting this to use the DRMAA API to interface with an SGE
#       cluster. From what the docs describe, DRMAA should let us interface with
#       several different distributed system implementations (e.g. TORQUE, PBS,
#       SGE). This will probably simplify the code below, and make things more
#       portable.

class SgeJobScheduler(AbstractJobScheduler):
    """
    Wrapper around the Sun Grid Engine (SGE) scheduler. Provides methods
    for monitoring job status, as well as submitting and killing jobs.

    The following code was developed and tested on SGE v6.2u5.
    """

    # The max number of seconds to wait for a completed job to get valid output
    # from the qacct command, before declaring the job failed.
    _MAX_WAIT_FOR_QACCT_AFTER_JOB_COMPLETE = 10

    # Regular expression for matching qstat output format, including header.
    _SGE_QSTAT_OUTPUT_PATTERN = re.compile(r'''\s*job-ID\s+prior\s+name\s+user\s+state\s+submit/start at\s+queue\s+slots\s+ja-task-ID\s+\n.*''')

    # Regular expression for parsing individual qstat output lines
    _SGE_QSTAT_OUTPUT_LINE_PATTERN = re.compile(r'''\s*(?P<job_id>\d+?)\s+\S+\s+\S+\s+\S+\s+(?P<job_status>\S+?)\s+.*''')

    # Regular expression for recognizing and extracting SGE job IDs following qsub.
    _SGE_QSUB_OUTPUT_PATTERN = re.compile(r'Your job (?P<job_id>\d+?) \(".*"\) has been submitted')

    # Regular expressions for identifying successful job deletion. Using separate
    # regex patterns since they are all capturing the same variable name. At the
    # moment, this is not needed as the job ID in the qdel output is not used for
    # anything. This is something to simplify in the future.
    #   Returned when qdel run on a running job
    _SGE_QDEL_OUTPUT_PATTERN_1 = re.compile(r'.* has registered the job (?P<job_id>\d+?) for deletion')
    #   Returned when qdel run on a pending job
    _SGE_QDEL_OUTPUT_PATTERN_2 = re.compile(r'.* has deleted job (?P<job_id>\d+?)')
    #   Returned when qdel run on job a second time, but it has not been deleted yet.
    _SGE_QDEL_OUTPUT_PATTERN_3 = re.compile(r'job (?P<job_id>\d+?) is already in deletion')

    # Default command used to check status of all jobs.
    _DEFAULT_QSTAT_COMMAND = ('qstat {qstat_args}')

    # Default command used to submit job.
    _DEFAULT_QSUB_COMMAND = ('qsub -N \"{job_name}\"'
                             ' -V -cwd'
                             ' -pe smp {num_processors}'
                             ' -l h_vmem={mem_usage_in_mb}M')

    # Default command used to kill job.
    _DEFAULT_QDEL_COMMAND = ('qdel {qdel_args} {job_id}')

    # Default command used to get detailed accounting info for a given job.
    _DEFAULT_QACCT_COMMAND = ('qacct {qacct_args} -j {job_id}')

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

        jobid_to_run_state = SgeJobScheduler._run_and_parse_qstat()
        # Convert job_id to string, since keys in jobid_to_run_state are strings
        # also (and integer keys won't return anything).
        sge_job_status = jobid_to_run_state.get(str(job_id))

        if sge_job_status:
            '''
            Excerpt from the qstat man page above job stats in SGE (version GE 6.2u5):

            -d(eletion) indicates that a qdel(1) has been used to initiate job deletion.
            -t(ransfering) and r(unning) indicate that a job is about to be executed or
             is already executing
            -s(uspended), S(uspended) and T(hreshold) show that an already running jobs
             has been suspended. The s(uspended) state is caused by suspending the job via
             the qmod(1) command, the S(uspended) state indicates that the queue containing
             the job is suspended and therefore the job is also suspended and the T(hreshold)
             state shows that at least one suspend threshold of the corresponding queue was
             exceeded (see queue_conf(5)) and that the job has been suspended as a consequence.
            -R(estarted) indicates that the job was restarted. This can be caused by a job
             migration or because of one of the reasons described in the -r section of the
             qsub(1) command.
            -w(aiting) and h(old) only appear for pending jobs. The h(old) state indicates
             that a job currently is not eligible for execution due to a hold state assigned
             to it via qhold(1), qalter(1) or the qsub(1) -h option or that the job is waiting
             for completion of the jobs to which job dependencies have been assigned to the
             job via the -hold_jid or -hold_jid-ad options of qsub(1) or qalter(1).
            -E(rror) appears for pending jobs that couldn't be started due to job properties.
             The reason for the job error is shown by the qstat(1) -j job_list option.
            '''
            if sge_job_status.startswith("E") or sge_job_status.startswith("d") or \
               sge_job_status == "s" or sge_job_status == "S":
                job_status = "FAILED"
                # TODO: In SGE, I think jobs with these states can linger in the queue.
                #       Might be worth adding a qdel command to remove jobs that have
                #       entered this failure state.
            elif sge_job_status == "r" or sge_job_status == "t":
                job_status = "RUNNING"
            elif sge_job_status.endswith("w") or sge_job_status.startswith("h"):
                job_status = "PENDING"

        else:

            #Check list of completed jobs
            jobid_to_run_state = SgeJobScheduler._run_and_parse_qstat(additional_args="-s z")

            if jobid_to_run_state.get(str(job_id)):
                #Use qacct to find the exit code for this job.
                job_exit_code = SgeJobScheduler._get_exit_code_from_qacct(job_id)

                '''
                There is often a short lag period bewteen when a job is finished
                (i.e. vieable in the list of completed jobs) and when the job
                accounting info is available (output of qacct). So if the job
                status is checked during this window, it's possible a job that
                has completed correctly will return an error status, because
                output from the qacct command isn't available yet. The current
                hack is to wait for some time interval if there is no qacct outpout,
                re-check the qacct output, and then proceed.

                NOTE: The only situations I've seen so far where job accounting
                      is never reported for a completed job is if the job was
                      killed while it was in the pending state. It seems the job
                      must at least get to the running state before it will generate
                      any accounting information.
                '''
                if job_exit_code == -1:
                    time.sleep(SgeJobScheduler._MAX_WAIT_FOR_QACCT_AFTER_JOB_COMPLETE)
                    job_exit_code = SgeJobScheduler._get_exit_code_from_qacct(job_id)

                if job_exit_code == 0:
                    job_status = "COMPLETED"
                elif job_exit_code > 0:
                    job_status = "FAILED"

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
            using the "-o" argument. Default: None.
        stderr_logfile : string
            Full path to file where job stderr should be stored. Specified using
            the "-e" argument. Default: None.
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

        qsub_command = SgeJobScheduler._DEFAULT_QSUB_COMMAND.format(job_name=job_name,
                                                                    num_processors=num_processors,
                                                                    mem_usage_in_mb=memory_in_mb)
        if stdout_logfile:
            qsub_command += f" -o {stdout_logfile}"
        if stderr_logfile:
            qsub_command += f" -e {stderr_logfile}"

        qsub_result = subprocess.run(' '.join([qsub_command, additional_args, job_command]),
                                     shell=True, check=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, encoding="ascii")
        if SgeJobScheduler._SGE_QSUB_OUTPUT_PATTERN.match(qsub_result.stdout):
            job_id = SgeJobScheduler._SGE_QSUB_OUTPUT_PATTERN.match(qsub_result.stdout).group('job_id')

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

        qdel_command = SgeJobScheduler._DEFAULT_QDEL_COMMAND.format(job_id=job_id, qdel_args=additional_args)
        # Note, set check=False here, since qdel command exits with an error code
        # if the given job does not exist or has already completed. These cases
        # are handled by the code afterward, so we don't need to check here (which
        # causes the whole run() statement to throw an exception).
        qdel_result = subprocess.run(qdel_command, shell=True, check=False,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     encoding="ascii")
        qdel_message = qdel_result.stdout.rstrip()

        if SgeJobScheduler._SGE_QDEL_OUTPUT_PATTERN_1.match(qdel_message) or \
           SgeJobScheduler._SGE_QDEL_OUTPUT_PATTERN_2.match(qdel_message) or \
           SgeJobScheduler._SGE_QDEL_OUTPUT_PATTERN_3.match(qdel_message):
            qdel_status = True

        #TODO: Wait and check job status to make sure job was really killed. Could
        #      escalate to different qdel arguments if the job fails to die, or
        #      simply throw an exception and let the calling code decided how
        #      to handle this.
        return qdel_status

    @staticmethod
    def _run_and_parse_qstat(additional_args=""):
        """
        Helper function that uses the qstat command to get the run status for all
        jobs in the SGE queue. Returns a dictionary of job ID's mapped to their
        run state codes. This is necessary because qstat does not return a job's
        run state when it is run on a single job (i.e. 'qstat -j <job_id>').

        Note: if a job has finished running, it will be missing from the default
        qstat command. To get a list of completed jobs, run this function with
        additional_args='-s z'. The states listed for the completed jobs *DO NOT*
        reflect whether the jobs completed successfully or with errors.

        Parameters
        ----------
        additional_args : string
            Additional arguments to provide to the qstat command. Settin this to
            '-s z' will return the list of completed jobs. Default: empty string.

        Returns
        -------
        dict
            Dictionary of job ID's mapped to their run state codes.

        """
        jobid_to_run_state = {}

        qstat_command = SgeJobScheduler._DEFAULT_QSTAT_COMMAND.format(qstat_args=additional_args)
        qstat_result = subprocess.run(qstat_command, shell=True, check=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      encoding="ascii")

        if SgeJobScheduler._SGE_QSTAT_OUTPUT_PATTERN.match(qstat_result.stdout):
            for line in qstat_result.stdout.splitlines():
                #This will skip header line, since "job-ID" column label doesn't
                #consist entire of digits.
                if SgeJobScheduler._SGE_QSTAT_OUTPUT_LINE_PATTERN.match(line):
                    job_id = SgeJobScheduler._SGE_QSTAT_OUTPUT_LINE_PATTERN.match(line).group('job_id')
                    job_status = SgeJobScheduler._SGE_QSTAT_OUTPUT_LINE_PATTERN.match(line).group('job_status')
                    jobid_to_run_state[job_id] = job_status

        return jobid_to_run_state

    @staticmethod
    def _get_exit_code_from_qacct(job_id, additional_args=""):
        """
        Helper method that runs the qacct command for the given job ID, and parses
        the output to find and return the job's exit code. This is necessary because
        the qstat output does not indicate the state in which jobs completed. Also,
        the job-specific qstat output (which contains more detailed information) is
        only available for running jobs.

        Note: qacct output is only available for completed commands.

        Parameters
        ----------
        job_id : string
            Unique SGE job id.
        additional_args : string
            Additional arguments to provide to the qacct command. Default: empty string.

        Returns
        -------
        string
            Exit code for the job, as reported by the qacct command. If the job
            was not found by the qacct command, returns "-1".

        """
        exit_code = -1

        qacct_command = SgeJobScheduler._DEFAULT_QACCT_COMMAND.format(job_id=job_id,
                                                                      qacct_args=additional_args)
        # Note, set check=False here, since qacct command exits with an error code
        # if accoungint for the given job does not exist. This case is handled by
        # the code afterward, so we don't need to check here (which causes the whole
        # run() statement to throw an exception).
        qacct_results = subprocess.run(qacct_command, shell=True, check=False,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       encoding="ascii")

        if qacct_results.stdout != f'error: job id {job_id} not found':
            for line in qacct_results.stdout.splitlines():
                if line.startswith("exit_status"):
                    key, exit_code = line.split(None, 1)
                    exit_code = int(exit_code)

        return exit_code
