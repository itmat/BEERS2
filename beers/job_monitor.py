import os
import sys
import time
import re
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

    def __init__(self, output_directory_path, process_list, samples, dispatcher_mode, step_name):
        """
        Initialize the monitor to track a specific set of jobs/processes running on
        a list of corresponding samples.

        Parameters
        ----------
        output_directory_path : string
            Path to data directory where job/process output is being stored.
        process_list : dict
            Mapping of job/process IDs to their associated sample IDs. Job
            monitoring involves tracking which jobs/processes are still running
            (requires job/process IDs), as well as parsing log files (requires
            sample IDs).
        samples : dict
            Sample information (IDs, input files). This might not be necessary if
            the only info I need is the sample ID.
        dispatcher_mode : string
            The mode used to submit the jobs/processes. Currently supports
            {",".join(SUPPORTED_DISPATCHER_MODES)}.
        step_name : string
            Name of the step in the pipeline being monitored. Ideally this code
            should be agnostic to the step, but there are some steps that will
            have different monitoring requirements (like STAR alignement, vs
            monitoring a molecule packet making its way through the library_prep
            pipeline). I might eventually be able to do away with this, or
            generalize the code further.
        """

        #TODO: add checks to make sure process_list and samples lengths are
        #      both equal and non-zero.
        self.output_directory = output_directory_path
        self.log_directory = os.path.join(self.output_directory, CONSTANTS.LOG_DIRECTORY_NAME)
        self.data_directory = os.path.join(self.output_directory, CONSTANTS.DATA_DIRECTORY_NAME)
        self.process_list = process_list
        self.samples = samples
        self.step_name = step_name
        #Tracks which of the submitted jobs have stalled or failed and ultimately
        #will require resubmission. I might want to separate stalled jobs from
        #failed jobs, or maybe I just include flags in the dictionary which
        #indicate why the job is listed as failed.
        self.resubmission_list = {}
        self.completed_list = {}

        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise BeersException(f'{dispatcher_mode} is not a supported mode.\n'
                                 'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.dispatcher_mode = dispatcher_mode


    def is_processing_complete(self):
        #Check run status of each job/process. If both the process list and the
        #resubmission lists are empty, return true. Otherwise return false.
        #Note, I need to force python to create a copy of the process_list
        #so that if/when the code below removes processes from the process_lists
        #it won't cause python to throw a "dictionary changed size during iteration"
        #error.
        for process_id, sample_id in dict(self.process_list).items():
            self.check_job_status(process_id, sample_id)
        return True if len(self.process_list) == 0  and len(self.resubmission_list) == 0 else False

        #TODO: Solve problem with modification of process_list while iterating.
        #      If the only modifications are happening within the job_monitor
        #      class, this should be ok. However, if outside processes can modify
        #      the process_list (by calling resubmit_process, for example), this
        #      could cause problems if the modification happens while another
        #      piece of code is iterating over one of the lists. I could try to
        #      solve this by using the various process lists to buffer changes
        #      to the process list and those modifications are only carried out
        #      within the job_monitor class (e.g. after the is_processing_complete
        #      function finishes).

    def mark_process_completed(self, process_id, sample_id):
        self.completed_list[process_id] = sample_id
        del self.process_list[process_id]
        print("Job's done and added to completed list")

    def mark_process_for_resubmission(self, process_id, sample_id):
        self.resubmission_list[process_id] = sample_id
        #Remove job from the process list so it can be replaced with a new entry
        #in the process list rollowing resubmission.
        del self.process_list[process_id]

    def resubmit_process(self, process_id, sample_id):
        #TODO: Add check to make sure resubmitted process isn't already in the list
        #      in the process_list. Probably best to also verify the process ID
        #      is also in the resubmission list.
        self.process_list[process_id] = sample_id
        del self.resubmission_list[process_id]

    def check_job_status(self, process_id, sample_id):

        if self.dispatcher_mode == "serial" or self.dispatcher_mode == "parallel":
            pass
        elif self.dispatcher_mode == "lsf":

            result = subprocess.run(' '.join([f"bjobs {process_id}"]), shell=True, check=True, stdout=subprocess.PIPE, encoding="ascii")

            #Here's some code just using string split to try to get job status
            #Skip first line bjobs output, since it just contains the header info.
            #job_status = result.stdout.split("\n")[1].split()[2]

            if Monitor.lsf_bjobs_output_pattern.match(result.stdout):
                lsf_job_status = Monitor.lsf_bjobs_output_pattern.match(result.stdout).group("job_status")


                #Job still waiting or is currently running, so we don't care.
                #TODO: If jobs remain in either state for too long, resubmit,
                #      or check output/log files (in the case of RUN) for their
                #      their last update. If too much time has passed during
                #      an update, might need to resubmit jobs.
                if lsf_job_status == "PEND" or lsf_job_status == "RUN" or lsf_job_status == "WAIT":
                    pass
                elif lsf_job_status == "EXIT":
                    self.mark_process_for_resubmission(process_id, sample_id)
                elif lsf_job_status == "DONE":
                    #Check output files
                    if self.step_name == "GenomeAlignmentStep":
                        aligner_log_file_path = os.path.join(self.data_directory, f"sample{sample_id}", "genome_alignment.Log.progress.out")
                        #Read last line in aligner log file
                        with open(aligner_log_file_path, "r") as aligner_log_file:
                            line = ""
                            for line in aligner_log_file:
                                line = line.rstrip()
                            if line == "ALL DONE!":
                                self.mark_process_completed(process_id, sample_id)
                            else:
                                self.mark_process_for_resubmission(process_id, sample_id)
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
                self.mark_process_for_resubmission(process_id, sample_id)
        else:
            #TODO: generalize code, implement code for other steps, or both.
            raise NotImplementedError()
