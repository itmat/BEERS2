import os
import sys
import time
import re
import subprocess
from beers.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES
from beers.beers_exception import BeersException

class Job:
    """
    Wrapper around subprocesses executed throughout the pipeline. Contains
    methods for determining the job's run status and lists dependencies on any
    other jobs.
    """

    #Regular expression for parsion bjobs output (including header)
    lsf_bjobs_output_pattern = re.compile(r'''JOBID\s+USER\s+STAT\s+QUEUE\s+FROM_HOST\s+EXEC_HOST\s+JOB_NAME\s+SUBMIT_TIME\n(?P<job_id>\d+?)\s+\S+\s+(?P<job_status>\S+?)\s+.*''')

    def __init__(self, job_id, sample_id, step_name, output_directory_path,
                 dispatcher_mode, system_id=None):
        """
        Initialize job to track the status of a step/operation running on the
        current sample.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this specific job.
        sample_id : string
            ID of sample object (created by Controller class) associated with
            job.
        step_name : string
            Name of the step in the pipeline being monitored. Ideally this code
            should be agnostic to the step, but there are some steps that will
            have different monitoring requirements (like STAR alignement, vs
            monitoring a molecule packet making its way through the library_prep
            pipeline). I might eventually be able to do away with this, or
            generalize the code further.
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        dispatcher_mode : string
            The mode used to submit the jobs/processes. Currently supports
            {",".join(SUPPORTED_DISPATCHER_MODES)}.
        system_id : string
            System-level identifier for the running job. For example, if this
            job was submitted to a job schedule  like LSF or SGE, this will be
            the external identifier assigned by the LSF/SGE system. If the job
            has not been submitted to the system yet (e.g. it is waiting for one
            of its dependencies to finish), then the system_id should be set to
            "None" [default].
        """

        self.job_id = job_id
        self.sample_id = sample_id
        self.system_id = system_id
        self.step_name = step_name
        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.data_directory = os.path.join(self.output_directory, CONSTANTS.DATA_DIRECTORY_NAME)

        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise BeersException(f'{dispatcher_mode} is not a supported mode.\n'
                                 'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.dispatcher_mode = dispatcher_mode

        self.dependency_list = []

    def add_dependency(self, dependency_job_id):
        """
        Add given id to job's dependency list, if it is not already there.

        Parameters
        ----------
        dependency_job_id : string
            Internal BEERS ID to add to dependency list.

        """
        if dependency_job_id not in self.dependency_list:
            self.dependency_list.append(dependency_job_id)



    def check_job_status(self):
        """
        Determine job's current run status based on system's job handler status
        and the job's output files.

        Returns
        -------
        string
            SUBMITTED - job submitted to system (might be running or waiting).
            FAILED - job finished with error status or incomplete output files.
            STALLED - job running without any change in output files for longer than threshold time.
            COMPLETED - job finished successfully with complete output files.
            WAITING_FOR_DEPENDENCY - job not submitted and waiting for dependency to complete.
        """

        job_status = "SUBMITTED"

        if self.system_id is None:
            job_status = "WAITING_FOR_DEPENDENCY"
        elif self.dispatcher_mode == "serial" or self.dispatcher_mode == "parallel":
            pass
        elif self.dispatcher_mode == "lsf":

            result = subprocess.run(' '.join([f"bjobs {self.system_id}"]), shell=True, check=True, stdout=subprocess.PIPE, encoding="ascii")

            #Here's some code just using string split to try to get job status
            #Skip first line bjobs output, since it just contains the header info.
            #job_status = result.stdout.split("\n")[1].split()[2]

            if Job.lsf_bjobs_output_pattern.match(result.stdout):
                lsf_job_status = Job.lsf_bjobs_output_pattern.match(result.stdout).group("job_status")


                #Job still waiting or is currently running, so we don't care.
                #TODO: If jobs remain in either state for too long, resubmit,
                #      or check output/log files (in the case of RUN) for their
                #      their last update. If too much time has passed during
                #      an update, might need to resubmit jobs.
                if lsf_job_status == "PEND" or lsf_job_status == "RUN" or lsf_job_status == "WAIT":
                    pass
                elif lsf_job_status == "EXIT":
                    job_status = "FAILED"
                elif lsf_job_status == "DONE":
                    #Check output files
                    if self.step_name == "GenomeAlignmentStep":
                        aligner_log_file_path = os.path.join(self.data_directory, f"sample{self.sample_id}", "genome_alignment.Log.progress.out")
                        #Read last line in aligner log file
                        with open(aligner_log_file_path, "r") as aligner_log_file:
                            line = ""
                            for line in aligner_log_file:
                                line = line.rstrip()
                            if line == "ALL DONE!":
                                job_status = "COMPLETED"
                            else:
                                job_status = "FAILED"
                else:
                    #TODO: Handle all other possible status messages from bjobs. Search
                    #      the bjobs manpage for "JOB STATUS" to find the full list and
                    #      explanation of all job status values.
                    raise NotImplementedError()

            #TODO: Right now I'm just assuming anything that isn't found by bjobs
            #      is a failed job and should be restarted. Strictly speaking,
            #      I'm not sure if this is a fair assumption to make, especially
            #      if something caused the LSF system's job histroy to reset. This
            #      should probably also check some output and log files to verify
            #      that something did go wrong.
            else:
                job_status = "FAILED"
        else:
            #TODO: generalize code, implement code for other steps, or both.
            raise NotImplementedError()

        return job_status
