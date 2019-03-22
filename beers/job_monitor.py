import os
import re
import importlib
import subprocess
from beers.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES
from beers.beers_exception import BeersException

class Monitor:
    """
    The class monitors the status of various subprocesses running throughout the
    pipeline. It checks for jobs that are pending, running, stalled, or halted
    (either due to success or error/failure).
    """

    #Regular expression for parsion bjobs output (including header)
    lsf_bjobs_output_pattern = re.compile(r'''JOBID\s+USER\s+STAT\s+QUEUE\s+FROM_HOST\s+EXEC_HOST\s+JOB_NAME\s+SUBMIT_TIME\n(?P<job_id>\d+?)\s+\S+\s+(?P<job_status>\S+?)\s+.*''')

    #Regex for recognizing and extracting lsf job IDs following submission.
    lsf_bsub_output_pattern = re.compile(r'Job <(?P<job_id>\d+?)> is submitted .*')

    def __init__(self, output_directory_path, dispatcher_mode):
        """
        Initialize the monitor to track a specific set of jobs/processes running on
        a list of corresponding samples.

        Parameters
        ----------
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        dispatcher_mode : string
            The mode used to submit the jobs/processes. Currently supports
            {",".join(SUPPORTED_DISPATCHER_MODES)}.
        """

        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.pending_list = {}
        self.running_list = {}
        #Tracks which of the submitted jobs have stalled or failed and ultimately
        #will require resubmission. I might want to separate stalled jobs from
        #failed jobs, or maybe I just include flags in the dictionary which
        #indicate why the job is listed as failed.
        self.resubmission_list = {}
        self.completed_list = {}

        #Stores list of samples in dictionary indexed by sample ID.
        self.samples_by_ids = {}

        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise BeersException(f'{dispatcher_mode} is not a supported mode.\n'
                                 'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.dispatcher_mode = dispatcher_mode


    def is_processing_complete(self):
        #Check run status of each job. If the running job list, the pending list,
        #and the resubmission lists are empty, return true. Otherwise return false.
        #Note, I need to force python to create a copy of the running_list so
        #that if/when the code below removes bjos from the running_list it won't
        #cause python to throw a "dictionary changed size during iteration" error.
        for job_id, job in dict(self.running_list).items():
            job_status = job.check_job_status()
            if job_status == "FAILED":
                self.mark_job_for_resubmission(job_id)
            elif job_status == "COMPLETED":
                self.mark_job_completed(job_id)
        #TODO: Check pending_list's dependencies to see if they need to move to running_list.
        return True if (not self.running_list and
                        not self.pending_list and
                        not self.resubmission_list) else False

        #TODO: Solve problem with modification of running_list while iterating.
        #      If the only modifications are happening within the job_monitor
        #      class, this should be ok. However, if outside processes can modify
        #      the running_list (by calling resubmit_process, for example), this
        #      could cause problems if the modification happens while another
        #      piece of code is iterating over one of the lists. I could try to
        #      solve this by using the various process lists to buffer changes
        #      to the process list and those modifications are only carried out
        #      within the job_monitor class (e.g. after the is_processing_complete
        #      function finishes).

    def mark_job_completed(self, job_id):
        self.completed_list[job_id] = self.running_list[job_id]
        del self.running_list[job_id]

    def mark_job_for_resubmission(self, job_id):
        self.resubmission_list[job_id] = self.running_list[job_id]
        #Remove job from the process list so it can be replaced with a new entry
        #in the process list rollowing resubmission.
        del self.running_list[job_id]

    def submit_new_job(self, job_id, job_command, sample, step_name, validation_attributes,
                       output_directory_path, system_id=None, dependency_list=None):
        """
        Create a Job given the list of job attributes and add it to the running
        list if it has a system_id. If it has no system_id or has dependencies,
        add job to the pending list.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this job.
        job_command : string
            Full command to execute job when run from the command line. This will
            likely be the output of the StepName.get_commandline_call() function.
        sample : Sample
            Sample object associated with the job. Will be added to the dictionary
            of samples stored in the Monitor if it's not already there.
        step_name : string
            Name of the step in the pipeline associated with monitored job.
        validation_attributes : dict
            Dictionary of attribute names and values required to validate the
            output of this specific job. This will most often be generated by
            the get_validation_attributes() method of the corresponding step,
            and used by the is_output_valid() method to find and test the integrity
            of any output files or parameters.
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        system_id : string
            System-level identifier for the running job. For example, if this
            job was submitted to a job schedule  like LSF or SGE, this will be
            the external identifier assigned by the LSF/SGE system. If the job
            has not been submitted to the system yet (e.g. it is waiting for one
            of its dependencies to finish), then the system_id should be set to
            "None" [default].
        dependency_list : list
            List of BEERS job IDs the submitted job is dependent upon (i.e. this
            job will wait until all those on the dependency list have completed).
            If the job has no dependencies, this should be "None" or empty.
        """
        submitted_job = Job(job_id, job_command, sample.sample_id, step_name,
                            validation_attributes, output_directory_path,
                            self.dispatcher_mode, system_id, dependency_list)
        #TODO: Condsider whether to add a check to see if both a system_id and
        #      a dependency list is provided (and throw an expection if both or
        #      neither are). As it's coded here, if a job has a system_id, it
        #      will be added to the running_list, regardless of whether or not
        #      its dependencies are met (if it has any). This puts the onus on
        #      whichever code is submitting jobs to make sure the system_id is
        #      absent from jobs that should start in the pending queue. Relaxing
        #      any sort of test here does give other code greater flexibility in
        #      assigning jobs to the various job lists, which could be particularly
        #      useful when restoring the monitor following a system crash or
        #      restart (e.g. jobs in a dependency list may not have been placed
        #      in the completed list if they were from an earlier step that didn't
        #      require restarting).

        if job_id in self.running_list or job_id in self.pending_list:
            raise MonitorException(f'Submitted job is already in the list of running or pending\n',
                                   'jobs. To move a job from the pending to the running list, use\n',
                                   'the submit_pending_job() function\n')
        else:
            if system_id is not None:
                self.running_list[job_id] = submitted_job
            else:
                self.pending_list[job_id] = submitted_job

            if not sample.sample_id in self.samples_by_ids:
                self.samples_by_ids[sample.sample_id] = sample


    def submit_pending_job(self, job_id, new_system_id):
        """
        Update job with new system ID and move it from the pending list to
        the list of running jobs.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies the job to submit. This
            job ID must be present in the pending list.
        new_system_id : string
            New system-level ID assigned to the job during submission.
        """
        if not job_id in self.pending_list:
            raise MonitorException(f'Job missing from the list of pending jobs.\n')
        elif job_id in self.running_list or job_id in self.resubmission_list:
            raise MonitorException(f'Job is already in the list of running jobs or ',
                                   'jobs marked for resubmission.\n')
        else:
            job = self.pending_list[job_id]
            job.system_id = new_system_id
            self.running_list[job_id] = job
            del self.pending_list[job_id]

    def resubmit_job(self, job_id, new_system_id):
        """
        Update job with new system ID and move it from the resubmission list to
        the list of running jobs. Also increment job's resubmission counter.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies the job to resubmit. This
            job ID must be present in the resubmission list.
        new_system_id : string
            New system-level ID assigned to the job during resubmission.
        """
        if not job_id in self.resubmission_list:
            raise MonitorException(f'Resubmitted job missing from the list of jobs ',
                                   'marked for resubmission.\n')
        elif job_id in self.running_list or job_id in self.pending_list:
            raise MonitorException(f'Resubmitted job is already in the list of ',
                                   'running or pending jobs.\n')
        else:
            job = self.resubmission_list[job_id]
            job.system_id = new_system_id
            job.resubmission_counter += 1
            self.running_list[job_id] = job
            del self.resubmission_list[job_id]

    def get_sample(self, sample_id):
        """
        Helper function returns sample object given a sample_id.
        """
        return self.samples_by_ids[sample_id]

    def are_dependencies_satisfied(self, job_id):
        """
        Given a job_id, check to see if all of its dependencies are satisfied
        (i.e. all dependencies are in the completed list).

        Parameters
        ----------
        job_id : string
            Internal BEERS ID of a job present in the pending list.

        Returns
        -------
        boolean
            True - all of the given job's depenencies are in the completed list
                   or if the job has no dependencies.
            False - the job has dependnecies that are not in the completed list.
        """
        #TODO: If we need to optimize in the future, we could make the dependency
        #      list a set, rather than a list, which would save us the set conversion
        #      operation.
        job = self.pending_list[job_id]
        return set(job.dependency_list).issubset(self.completed_list.keys())


class MonitorException(BeersException):
    pass


class Job:
    """
    Wrapper around subprocesses executed throughout the pipeline. Contains
    methods for determining the job's run status and lists dependencies on any
    other jobs.
    """

    #Regular expression for parsion bjobs output (including header)
    LSF_BJOBS_OUTPUT_PATTERN = re.compile(r'''JOBID\s+USER\s+STAT\s+QUEUE\s+FROM_HOST\s+EXEC_HOST\s+JOB_NAME\s+SUBMIT_TIME\n(?P<job_id>\d+?)\s+\S+\s+(?P<job_status>\S+?)\s+.*''')
    #List of valid outputs from job status-reporting methods. This is a guide
    #for future development and users who wish to add custom status-checking
    #methods.
    JOB_STATUS_OUTPUTS = ['SUBMITTED',              #job submitted to system (might be running or waiting).
                          'FAILED',                 #job finished with error status or incomplete output files.
                          'STALLED',                #job running without any change in output files for longer than threshold time.
                          'COMPLETED',              #job finished successfully with complete output files.
                          'WAITING_FOR_DEPENDENCY'] #job not submitted and waiting for dependency to complete.

    def __init__(self, job_id, job_command, sample_id, step_name, validation_attributes,
                 output_directory_path, dispatcher_mode, system_id=None,
                 dependency_list=None):
        """
        Initialize job to track the status of a step/operation running on the
        current sample.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this specific job.
        job_command : string
            Full command to execute job when run from the command line. This will
            likely be the output of the StepName.get_commandline_call() function.
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
        validation_attributes : dict
            Dictionary of attribute names and values required to validate the
            output of this specific job. This will most often be generated by
            the get_validation_attributes() method of the corresponding step,
            and used by the is_output_valid() method to find and test the integrity
            of any output files or parameters.
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
        dependency_list: list
            List of BEERS job IDs that this job is dependent upon (i.e. this
            job will wait until all those on the dependency list have completed).
            Empty list or "None" if there are no dependencies [default].
        """

        self.job_id = job_id
        self.sample_id = sample_id
        self.job_command = job_command
        self.step_name = step_name
        self.validation_attributes = validation_attributes
        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.data_directory = os.path.join(self.output_directory, CONSTANTS.DATA_DIRECTORY_NAME)

        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise BeersException(f'{dispatcher_mode} is not a supported mode.\n'
                                 'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.dispatcher_mode = dispatcher_mode

        self.system_id = system_id
        self.dependency_list = set()
        if dependency_list:
            self.add_dependencies(dependency_list)

        #Tracks number of times job has been resubmitted. Can use this to kill
        #entire pipeline if single jobs fail too many times (might indicate
        #something other than system instability is the cause for the failure
        #and should be remedied.
        self.resubmission_counter = 0


    def add_dependencies(self, dependency_job_ids):
        """
        Add given ids to job's dependency list, if they are not already there.

        Parameters
        ----------
        dependency_job_ids : list
            Internal BEERS IDs to add to dependency list.

        """
        for job_id in dependency_job_ids:
            self.dependency_list.add(job_id)

    def check_job_status(self):
        """
        Determine job's current run status based on system's job handler status
        and the job's output files.

        Returns
        -------
        string
            One of the following:
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
            job_status = "COMPLETED"
        elif self.dispatcher_mode == "lsf":

            result = subprocess.run(' '.join([f"bjobs {self.system_id}"]), shell=True, check=True, stdout=subprocess.PIPE, encoding="ascii")

            #Here's some code just using string split to try to get job status
            #Skip first line bjobs output, since it just contains the header info.
            #job_status = result.stdout.split("\n")[1].split()[2]

            if Job.LSF_BJOBS_OUTPUT_PATTERN.match(result.stdout):
                lsf_job_status = Job.LSF_BJOBS_OUTPUT_PATTERN.match(result.stdout).group("job_status")


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
                    #TODO: To clean up the code, we should probably create
                    #      separate methods to check the status of various steps,
                    #      rather than build them all in here. This code will
                    #      check the name of the step and call the appropriate
                    #      one of these methods.

                    #Check output files
                    if self.step_name == "GenomeAlignmentStep":
                        aligner_log_file_path = os.path.join(self.data_directory, f"sample{self.sample_id}", "genome_alignment.Log.progress.out")
                        #Read last line in aligner log file
                        line = ""
                        with open(aligner_log_file_path, "r") as aligner_log_file:
                            for line in aligner_log_file:
                                line = line.rstrip()
                        if line == "ALL DONE!":
                            job_status = "COMPLETED"
                        else:
                            job_status = "FAILED"
                    elif self.step_name == "GenomeBamIndexStep":
                        index_file_path = os.path.join(self.data_directory, f"sample{self.sample_id}", "genome_alignment.Aligned.sortedByCoord.out.bam.bai")
                        if os.path.isfile(index_file_path):
                            job_status = "COMPLETED"
                        else:
                            job_status = "FAILED"
                    elif self.step_name == "VariantsFinderStep":
                        module = importlib.import_module(f'.variants_finder', package="beers.expression")
                        VariantsFinderStep = getattr(module, self.step_name)
                        if VariantsFinderStep.is_output_valid(self.validation_attributes):
                            job_status = "COMPLETED"
                        else:
                            job_status = "FAILED"
                    elif self.step_name == "IntronQuantificationStep":
                        # TODO: do proper validation here
                        from beers.expression.intron_quant import IntronQuantificationStep
                        status = IntronQuantificationStep.is_output_valid({})
                        if status:
                            job_status = "COMPLETED"
                        else:
                            job_status = "FAILED"
                    else:
                        raise NotImplementedError()
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
