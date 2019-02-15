import os
import re
from beers.job import Job
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

    def submit_new_job(self, job_id, sample, step_name, output_directory_path,
                       system_id=None, dependency_list=None):
        """
        Create a Job given the list of job attributes and add it to the running
        list if it has a system_id. If it has no system_id or has dependencies,
        add job to the pending list.

        Parameters
        ----------
        job_id : string
            Internal BEERS ID that uniquely identifies this job.
        sample : Sample
            Sample object associated with the job. Will be added to the dictionary
            of samples stored in the Monitor if it's not already there.
        step_name : string
            Name of the step in the pipeline associated with monitored job.
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
        submitted_job = Job(job_id, sample.sample_id, step_name,
                            output_directory_path, self.dispatcher_mode,
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
