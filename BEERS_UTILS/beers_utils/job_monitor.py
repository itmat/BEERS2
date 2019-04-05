import os
import sys
import time
from beers_utils.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES
from beers_utils.general_utils import BeersUtilsException
import beers_utils.job_scheduler_provider

class JobMonitor:
    """
    The class monitors the status of various subprocesses running throughout the
    pipeline. It checks for jobs that are pending, running, stalled, or halted
    (either due to success or error/failure).
    """

    #Hash mapping step names (keys) to AbstractPipelineStep objects (values).
    #Used in code checking job output validity. Provides generalized way of
    #retrieving job-specific Class objects to run static methods when validating
    #job output.
    PIPELINE_STEPS = {}

    def __init__(self, output_directory_path, scheduler_name, max_resub_limit=3):
        """
        Initialize the monitor to track a specific set of jobs/processes running on
        a list of corresponding samples.

        Parameters
        ----------
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        scheduler_name : string
            The mode used to submit the jobs/processes to the scheduler.
        max_resub_limit : int
            The maximum number of times a job can be resubmitted before the
            pipeline halts.

        """
        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.max_resub_limit = max_resub_limit
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

        self.scheduler_name = scheduler_name
        self.job_scheduler = beers_utils.job_scheduler_provider.SCHEDULERS.get(scheduler_name)

    def is_processing_complete(self):
        """
        Check if all jobs have finished processing.

        Returns
        -------
        boolean
            True  - The running job list, pending list, and resubmission lists
                    are all empty.
            False - At least one job remains in the running job list, pending list,
                    or resubmission list.

        """
        # TODO: Could we merge this function with monitor_until_all_jobs_completed()?
        #       Would there every be any need to run is_processing_complete() alone?

        #Note, I need to force python to create a copy of the running_list so
        #that if/when the code below removes jobs from the running_list it won't
        #cause python to throw a "dictionary changed size during iteration" error.
        for job_id, job in dict(self.running_list).items():
            job_status = job.check_job_status(self.job_scheduler)
            if job_status == "FAILED":
                self.mark_job_for_resubmission(job_id)
            elif job_status == "COMPLETED":
                self.mark_job_completed(job_id)

        print(f"Running jobs:{len(self.running_list)} | "
              f"Pending jobs:{len(self.pending_list)} | "
              f"Resub jobs:{len(self.resubmission_list)} | "
              f"Completed jobs:{len(self.completed_list)}")

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
        """
        Move job from running queue to completed queue.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this job.

        """
        self.completed_list[job_id] = self.running_list[job_id]
        del self.running_list[job_id]

    def mark_job_for_resubmission(self, job_id):
        """
        Move job from running queue to resubmission queue.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this job.

        """
        self.resubmission_list[job_id] = self.running_list[job_id]
        #Remove job from the process list so it can be replaced with a new entry
        #in the process list rollowing resubmission.
        del self.running_list[job_id]

    def monitor_until_all_jobs_completed(self, queue_update_interval=10):
        """
        Monitor until all jobs in the pending, running, and resubmission queues
        have completed. This method performs the following job queue operations:
            1) Submit pending jobs as their dependencies are satisfied
            2) Mark and resubmit failed jobs.
            3) Move jobs to completed queue as they finish.

        Parameters
        ----------
        queue_update_interval : int
            Number of second to wait after checking and updating all jobs on the
            queues, before checking again.

        """
        while not self.is_processing_complete():
            #Check for jobs requiring resubmission
            resubmission_jobs = self.resubmission_list.copy()
            if resubmission_jobs:
                print(f"--Resubmitting {len(resubmission_jobs)} jobs that failed/stalled.")
                for resub_job_id in resubmission_jobs.keys():
                    self.resubmit_job(resub_job_id)

            #Check if pending jobs have satisfied their dependencies
            pending_jobs = self.pending_list.copy()
            if pending_jobs:
                print(f"--Check {len(pending_jobs)} pending jobs for satisfied dependencies:")
                for pend_job_id in pending_jobs.keys():
                    # TODO: Consider moving this check inside the submit_pending_job
                    #       method. The advantage of keeping this separate is we
                    #       can force the submission of a pending job, regardless
                    #       of the status of its dependencies (could be useful when
                    #       restarting a job monitoring queue following a crash).
                    if self.are_dependencies_satisfied(pend_job_id):
                        self.submit_pending_job(pend_job_id)
            time.sleep(queue_update_interval)

    def submit_new_job(self, job_id, job_command, sample, step_name, scheduler_arguments,
                       validation_attributes, output_directory_path, system_id=None,
                       dependency_list=None):
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
            of samples stored in the JobMonitor if it's not already there.
        step_name : string
            Name of the step in the pipeline associated with monitored job.
        scheduler_arguments : dict
            Dictionary of arguments passed to the job scheduler when submitting
            the job. These will be passed on to the submit_job() method of the
            job scheduler.
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
                            scheduler_arguments, validation_attributes,
                            output_directory_path, self.scheduler_name,
                            system_id, dependency_list)
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
            raise JobMonitorException(f"Submitted job is already in the list of running or pending\n"
                                      f"jobs. To move a job from the pending to the running list, use\n"
                                      f"the submit_pending_job() function\n")
        else:
            if system_id is not None:
                self.running_list[job_id] = submitted_job
            else:
                self.pending_list[job_id] = submitted_job

            if not sample.sample_id in self.samples_by_ids:
                self.samples_by_ids[sample.sample_id] = sample


    def submit_pending_job(self, job_id):
        """
        submit job through the job scheduler and move it from the pending list to
        the list of running jobs.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies the job to submit. This
            job ID must be present in the pending list.

        """
        # TODO: We could probably combine the resubmit_job() and submit_pending_job()
        #       into a single method that also takes a job_type argument [pending,resub]
        #       and then behaves accordingly. Or maybe it simply infers what to do with
        #       the job based on which queue it's located in. There's a ton of overlap
        #       between these two functions right now.
        if not job_id in self.pending_list:
            raise JobMonitorException(f"Job missing from the list of pending jobs.\n")
        elif job_id in self.running_list or job_id in self.resubmission_list:
            raise JobMonitorException(f"Job is already in the list of running jobs or "
                                      f"jobs marked for resubmission.\n")
        else:
            job = self.pending_list[job_id]

            pend_sample = self.get_sample(job.sample_id)

            print(f"\tSubmitting {job.step_name} command to {self.scheduler_name} "
                  f"for sample {pend_sample.sample_name}.")

            #Use unpacking to provide arguments for job submission
            new_system_id = self.job_scheduler.submit_job(job_command=job.job_command,
                                                          **job.scheduler_arguments)
            if new_system_id == "ERROR":
                print(f"Job submission failed for {job.step_name}:\n",
                      f"   Job sample: {pend_sample.sample_name}\n",
                      f"   Scheduler parameters: {job.scheduler_arguments}\n",
                      f"   Job command: {job.job_command}\n",
                      file=sys.stderr)
                raise JobMonitorException(f"Job submission failed for {job.step_name}. "
                                          f"See expression pipeline log file for full details.")

            print(f"\tFinished submitting {job.step_name} command to "
                  f"{self.scheduler_name} for sample {pend_sample.sample_name}.")

            job.system_id = new_system_id
            self.running_list[job_id] = job
            del self.pending_list[job_id]

    def resubmit_job(self, job_id):
        """
        Resubmit job through the job scheduler and move it from the resubmission
        list to the list of running jobs. Also increment job's resubmission counter.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies the job to resubmit. This
            job ID must be present in the resubmission list.

        """
        # TODO: We could probably combine the resubmit_job() and submit_pending_job()
        #       into a single method that also takes a job_type argument [pending,resub]
        #       and then behaves accordingly. Or maybe it simply infers what to do with
        #       the job based on which queue it's located in. There's a ton of overlap
        #       between these two functions right now.
        if not job_id in self.resubmission_list:
            raise JobMonitorException(f"Resubmitted job missing from the list of jobs "
                                      f"marked for resubmission.\n")
        elif job_id in self.running_list or job_id in self.pending_list:
            raise JobMonitorException(f"Resubmitted job is already in the list of "
                                      f"running or pending jobs.\n")
        else:
            job = self.resubmission_list[job_id]

            if job.resubmission_counter == self.max_resub_limit:
                raise JobMonitorException(f"The {job_id} exceeded the maximum resubmission"
                                          f"limit of {self.max_resub_limit}.\n")

            resub_sample = self.get_sample(job.sample_id)

            print(f"\tSubmitting {job.step_name} command to {self.scheduler_name} "
                  f"for sample {resub_sample.sample_name}.")

            #Use unpacking to provide arguments for job submission
            new_system_id = self.job_scheduler.submit_job(job_command=job.job_command,
                                                          **job.scheduler_arguments)
            if new_system_id == "ERROR":
                print(f"Job submission failed for {job.step_name}:\n",
                      f"   Job sample: {resub_sample.sample_name}\n",
                      f"   Scheduler parameters: {job.scheduler_arguments}\n",
                      f"   Job command: {job.job_command}\n",
                      file=sys.stderr)
                raise JobMonitorException(f"Job submission failed for {job.step_name}.")

            print(f"\tFinished submitting {job.step_name} command to "
                  f"{self.scheduler_name} for sample {resub_sample.sample_name}.")

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


class Job:
    """
    Wrapper around subprocesses executed throughout the pipeline. Contains
    methods for determining the job's run status, lists dependencies on any
    other jobs, and all other information required to submit/resubmit the
    job to the scheduler.
    """

    #List of valid outputs from job status-reporting methods. This is a guide
    #for future development and users who wish to add custom status-checking
    #methods.
    JOB_STATUS_OUTPUTS = ['SUBMITTED',              #job submitted to system (might be running or waiting).
                          'FAILED',                 #job finished with error status or incomplete output files.
                          'STALLED',                #job running without any change in output files for longer than threshold time.
                          'COMPLETED',              #job finished successfully with complete output files.
                          'WAITING_FOR_DEPENDENCY'] #job not submitted and waiting for dependency to complete.

    def __init__(self, job_id, job_command, sample_id, step_name, scheduler_arguments,
                 validation_attributes, output_directory_path, scheduler_name,
                 system_id=None, dependency_list=None):
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
        scheduler_arguments : dict
            Dictionary of arguments passed to the job scheduler when submitting
            the job. These will be passed on to the submit_job() function.
        validation_attributes : dict
            Dictionary of attribute names and values required to validate the
            output of this specific job. This will most often be generated by
            the get_validation_attributes() method of the corresponding step,
            and used by the is_output_valid() method to find and test the integrity
            of any output files or parameters.
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        scheduler_name : string
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
        self.scheduler_arguments = scheduler_arguments
        self.validation_attributes = validation_attributes
        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.data_directory = os.path.join(self.output_directory, CONSTANTS.DATA_DIRECTORY_NAME)

        if scheduler_name not in SUPPORTED_DISPATCHER_MODES:
            raise BeersUtilsException(f'{scheduler_name} is not a supported mode.\n'
                                      f'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.scheduler_name = scheduler_name

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

    def check_job_status(self, scheduler=None):
        """
        Determine job's current run status based on system's job handler status
        and the job's output files.

        Parameters
        ----------
        scheduler : AbstractJobScheduler
            Interface to the system's job scheduler currently tracking the job.
            If no job scheduler provided, the method assumes the job is running
            in locally in serial mode [default].

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
        elif self.scheduler_name == "serial" or self.scheduler_name == "parallel" or not scheduler:
            job_status = "COMPLETED"
        else:

            scheduler_job_status = scheduler.check_job_status(self.system_id)

            if scheduler_job_status == "RUNNING" or scheduler_job_status == "PENDING":
                job_status = "SUBMITTED"
            elif scheduler_job_status == "FAILED":
                job_status = "FAILED"
            elif scheduler_job_status == "COMPLETED":

                #TODO: To clean up the code, we should probably create
                #      separate methods to check the status of various steps,
                #      rather than build them all in here. This code will
                #      check the name of the step and call the appropriate
                #      one of these methods.

                #Check output files
                if self.step_name == "GenomeAlignmentStep" or \
                   self.step_name == "GenomeBamIndexStep" or \
                   self.step_name == "VariantsFinderStep" or \
                   self.step_name == "IntronQuantificationStep":

                    pipeline_step = JobMonitor.PIPELINE_STEPS[self.step_name]

                    if pipeline_step.is_output_valid(self.validation_attributes):
                        job_status = "COMPLETED"
                    else:
                        job_status = "FAILED"

                else:
                    raise NotImplementedError()

            #TODO: Right now I'm just assuming anything that isn't found by bjobs
            #      is a failed job and should be restarted. Strictly speaking,
            #      I'm not sure if this is a fair assumption to make, especially
            #      if something caused the LSF system's job histroy to reset. This
            #      should probably also check some output and log files to verify
            #      that something did go wrong.
            else:
                job_status = "FAILED"

        return job_status


class JobMonitorException(BeersUtilsException):
    pass
