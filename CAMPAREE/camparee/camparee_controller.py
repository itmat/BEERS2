import json
import time
import glob
import yaml
import numpy as np
import termcolor
import os
import sys
import traceback
import copy
import string
from datetime import datetime
from beers_utils.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES
from beers_utils.general_utils import GeneralUtils
from camparee.expression_pipeline import ExpressionPipeline,CampareeValidationException
from beers_utils.sample import Sample


class CampareeController:
    """
    The object essentially controls the flow of the pipeline.  The run_camparee.py command instantiates a controller and
    calls one of the controller methods depending upon the pipeline-stage requested.  The other methods in the class
    are helper methods.
    """

    """
    Get CAMPAREE root directory. Used to determin path the 'resources' subdirectory.
    """
    CAMPAREE_ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def __init__(self):
        """
        These attributes have no value at the time of instantiation.  They are populated by the helper methods.
        """
        self.controller_name = 'camparee_controller'
        self.resources_name = 'resources'
        self.controller_log_filename = 'camparee_controller.log'
        self.allowed_directory_name_set = set(string.ascii_letters + string.digits + "_")
        # The following attributes are defined following instantiation
        self.debug = False
        self.run_id = None
        self.dispatcher = None
        self.configuration = None
        self.resources = None
        self.seed = None
        self.output_directory_path = None
        self.input_samples = []

    def run_camparee_pipeline(self, args):
        """
        This is how run_camparee.py calls the camparee pipeline. This method reads
        the command line arguments, parses the config file,  and calls the necessary
        methods to run camparee.
        :param args: command line arguments
        """
        stage_name = "expression_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        if not self.validate_samples():
            raise CampareeValidationException("Sample data is not valid.  Please consult the standard error file"
                                              "for details.")
        #TODO: Once a dispatcher is more integrated with the epxpression pipeline,
        #we may want to move this check elsewhere.
        dispatcher_mode = ""
        if not args.dispatcher_mode:
            if not self.configuration['setup'].get('dispatcher_mode', None):
                raise CampareeValidationException('No dispatcher_mode given either on the command line'
                                                  ' or in the configuration file')
            dispatcher_mode = self.configuration['setup']['dispatcher_mode']
        else:
            dispatcher_mode = args.dispatcher_mode
        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise CampareeValidationException(f'{dispatcher_mode} is not a supported mode.\n'
                                              'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        self.assemble_input_samples()
        ExpressionPipeline.main(self.configuration, dispatcher_mode,
                                os.path.join(self.output_directory_path,stage_name),
                                self.input_samples)

    def perform_setup(self, args, stage_names):
        """
        This helper method sets up a number of attributes and behaviors in the controller.  Stacktraces are suppressed
        and only user friendly errors are shown when the debugger is off (just a command line arg right now).  The full
        configuration file data and run id are salted away and the random seed is set.  The initial output folder
        structure (excluding the subdirectory structure needed to accommodate large numbers of file) is created.  The
        output folder structure depends on the stage names.  Also, the controller log is started.
        :param args: The command line arguments
        :param stage_names: The stage names
        """
        self.debug = args.debug

        def exception_handler(exception_type, exception, tb):
            if args.debug:
                traceback.print_exception(exception_type, exception, tb)
            else:
                termcolor.cprint(f"ERROR: {exception}", 'red', file=sys.stderr)
        sys.excepthook = exception_handler
        self.retrieve_configuration(args.config)
        self.set_run_id(args.run_id)
        self.plant_seed(args.seed)
        self.create_output_folder_structure(stage_names)
        self.create_controller_log()
        os.sys.exit("Developer stop")

    def retrieve_configuration(self, configuration_file_path):
        """
        Helper method to parse the configuration file given by the path info a dictionary attached to the controller
        object.  For convenience, the portion of the configuration file that contains parametric data specific to the
        controller is set to a separate dictionary also attached to the controller.
        :param configuration_file_path: The absolute file path of the configuration file
        """
        with open(configuration_file_path, 'r+') as configuration_file:
            try:
                self.configuration = yaml.load(configuration_file, Loader=yaml.FullLoader)
            except yaml.YAMLError as ye:
                raise CampareeValidationException(ye)
        self.resources = self.configuration['resources']

    def set_run_id(self, run_id):
        """
        Helper method to add the run id to the controller for later use.  The run id, if any, on the command line
        takes precedence.  If no run id is present on the command line, the controller configuration data is searched
        for it.  If no run id is found in either place, an error is raised.
        :param run_id: The run id found on the command line, if any.
        """
        self.run_id = run_id
        if not self.run_id:
            if 'run_id' not in self.configuration['setup'] or not self.configuration['setup']['run_id']:
                raise CampareeValidationException('No run id was given either on the command line '
                                                  'or in the configuration file')
            self.run_id = self.configuration['setup']['run_id']
        if not set(self.run_id).issubset(self.allowed_directory_name_set):
            raise CampareeValidationException(f'The run id {self.run_id} contains disallowed characters.')


    def plant_seed(self, seed):
        """
        Helper method to add the seed to the controller for later use.  The seed, if any, on the command line, takes
        precedence.  If no seed is present on the command line, the controller configuration data is searched for it.
        If no seed is found there, one will be randomly generated.  The seed will be added to the controller log file
        created for this run so that the user may re-create the run exactly at a later date, assuming all else remains
        the same.
        :param seed: The seed value found on the command line, if any.
        """
        self.seed = seed
        if not self.seed:
            # If seed is not a config key or the seed is null, generate a seed internally
            if 'seed' not in self.configuration['setup'] or not self.configuration['setup']['seed']:
                self.seed = GeneralUtils.generate_seed()
            else:
                # Insure the config seed is an integer, otherwise throw an exception
                if not isinstance(self.configuration['setup']['seed'], int):
                    raise CampareeValidationException('The seed specified in the configuration file must be an integer')
                self.seed = self.configuration['setup']['seed']
        np.random.seed(self.seed)

    def create_output_folder_structure(self, stage_names):
        """
        Use the provided stage names, the run id and the top level output directory path from the configuration data
        to create the top level directory of a preliminary directory structure.  The attempt fails if either the top
        level output directory exists but either is not a directory or is a non-empty directory or if the user has
        insufficient permissions to create the directory. Created in the level directly below the top level output
        directory, are folders named after the stage names provided (i.e., controller, library_prep_pipeline,
        sequence_pipeline) and beneath each of these are data and log folders.  Additional subdirectories are created
        later to organize the numerous files exprected and avoid congestion.
        :param stage_names: names of folders directly below the top level output directory (e.g., controller,
        library_prep)
        """
        self.output_directory_path = os.path.join(self.configuration['output']['directory_path'],f"run_{self.run_id}")
        if not os.path.exists(self.output_directory_path):
            try:
                os.makedirs(self.output_directory_path, mode=0o0755, exist_ok=True)
            except PermissionError:
                raise CampareeValidationException(f"Insufficient permissions to create {self.output_directory_path}")
        else:
            if not os.path.isdir(self.output_directory_path):
                raise CampareeValidationException(f"The output directory path, {self.output_directory_path},"
                                                  f" is not a directory.")
            if os.listdir(self.output_directory_path):
                raise CampareeValidationException(f"The output directory path, {self.output_directory_path},"
                                                  f" must be empty.")
        for stage_name in stage_names:
            os.makedirs(os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME),
                        mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME),
                        mode=0o0755, exist_ok=True)

    def create_controller_log(self):
        """
        Generates a controller log containing the timestamp, seed, run id and current configuration data so that the
        user can replicate this run at a later date.
        """
        log_file_path = os.path.join(self.output_directory_path,
                                     self.controller_name,
                                     CONSTANTS.LOG_DIRECTORY_NAME,
                                     self.controller_log_filename)
        with open(log_file_path, 'w') as controller_log_file:
            timestamp = time.time()
            current_datetime = datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')
            controller_log_file.write(f"Run Id:\t{self.run_id}\n")
            controller_log_file.write(f"Run Timestamp (UTC):\t{current_datetime}\n")
            controller_log_file.write(f"Seed:\t{self.seed}\n")
            controller_log_file.write(f"Configuration:\n")
            json.dump(self.configuration, controller_log_file, indent=2)

    def validate_samples(self):
        """
        Iterates over the starting samples and verifies that their files can be found and that the gender designations,
        if any, are appropriate.
        :return: True is valid and false otherwise.
        """
        valid = True
        input_directory_path = self.configuration["input"]["directory_path"]
        self.input_samples = []
        for input_sample in self.configuration["input"]["data"].values():
            input_files = copy.copy(input_sample["fastq_files"])
            if "bam_file" in input_sample:
                input_files.append(input_sample["bam_file"])

            for filename in input_files:
                input_sample_file_path = os.path.join(input_directory_path, filename)
                if not os.path.exists(input_sample_file_path) or not os.path.isfile(input_sample_file_path):
                    print(f"The input sample file, {input_sample_file_path}, does not exist as a file.", file=sys.stderr)
                    valid = False
            gender = input_sample.get("gender", None)
            if gender:
                gender = gender.lower()
                if gender not in [CONSTANTS.MALE_GENDER, CONSTANTS.FEMALE_GENDER]:
                    print(f"The gender, {gender}, for input sample {input_sample['filename']} must be either"
                          f" {CONSTANTS.MALE_GENDER}, {CONSTANTS.FEMALE_GENDER} or not present.", file=sys.stderr)
                    valid = False
            else:
                print(f"The input sample, {input_sample['filenames']} has no gender specified.  Consequently, no"
                      f" gender specific chromosomes will be processed for this sample.")
        return valid



    def assemble_input_samples(self):
        """
        Creates a list of sample objects, attached to the controller, that represent those samples that are to be
        run in the expression pipeline.  If not running from the expression pipeline, this method is not used since
        the sample data is already contained in each packet.  For each sample, a unique combination of adapter sequences
        are provided.  The sample name is assumed to be that of the input filename without the extension.  Gender may
        or may not be provided in the configuration data.  If not set, the gender will be inferred by the expression
        pipeline.
        """
        input_directory_path = self.configuration["input"]["directory_path"]
        self.input_samples = []
        # TODO handle the situation where the adapter kit is not specified or not found
        # The kit is really only needed for library prep.  So if the expression pipeline does not generate
        # molecule packets, we could postpone this step until when that assembly occurs.  But we don't want to
        # make the addition to thousands of molecule packets after the fact.
        for sample_name, input_sample in self.configuration["input"]["data"].items():
            #sample_name = os.path.splitext(input_sample["filenames"][0])[0]
            fastq_file_paths = [os.path.join(input_directory_path, filename)
                                       for filename in input_sample["fastq_files"]]
            bam_file_path = os.path.join(input_directory_path, input_sample["bam_file"]) if "bam_file" in input_sample else ''
            gender = input_sample.get("gender", None)
            if gender:
                gender = gender.lower()
            self.input_samples.append(
                Sample(sample_id=Sample.next_sample_id,
                       sample_name=sample_name,
                       fastq_file_paths=fastq_file_paths,
                       adapter_sequences="",
                       bam_file_path=bam_file_path,
                       gender=gender))
            Sample.next_sample_id += 1
