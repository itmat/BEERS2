import json
import time
import pathlib
import numpy as np
import termcolor
import os
import shutil
import sys
import traceback
import copy
from datetime import datetime
from beers_utils.constants import CONSTANTS,SUPPORTED_SCHEDULER_MODES
from beers_utils.general_utils import GeneralUtils
from beers_utils.read_fasta import read_fasta
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.flowcell import Flowcell
from beers_utils.molecule_packet import MoleculePacket
from beers.dispatcher import Dispatcher
from beers.fast_q import FastQ
from beers.sam import SAM
from beers.auditor import Auditor
from beers.validator import Validator
import glob


class Controller:
    """
    The object essentially controls the flow of the pipeline.  The run_beers.py command instantiates a controller and
    calls one of the controller methods depending upon the pipeline-stage requested.  The other methods in the class
    are helper methods.  Additionally, the controller handles those tasks that cannot be distributed to isolated nodes.
    """

    def __init__(self):
        """
        These attributes have no value at the time of instantiation.  They are populated by the helper methods.
        """
        self.controller_name = 'controller'
        self.resources_name = 'resources'
        self.controller_log_filename = 'controller.log'
        # The following attributes are defined following instantiation
        self.debug = False
        self.run_id = None
        self.dispatcher = None
        self.configuration = None
        self.configuration_file_path = None
        self.resources = None
        self.controller_configuration = None
        self.flowcell = None
        self.seed = None
        self.output_directory_path = None
        self.input_samples = []

    def run_library_prep_pipeline(self, args):
        """
        This is how run_beers.py calls the library prep pipeline by itself.  This pipeline is complete and functional.
        Controller attributes are set up.  All molecule packet files are located.  The data and log output directories
        are identified and the data directory and number of packets to be processed are used to fully develop the
        output folder structure, which may have nested sub-directories in the case of a very large number of input
        molecule packets.  The step log directories are created.  The dispatcher is instantiated and finally, the
        dispatcher is run with the molecule_packet_file_paths provided.  Note that this pipeline stage assumes all
        molecule packets are immediately available.

        Note that the molecule packets output from this procedure have been pared down in size according to the
        flowcell retention percentage specified in the configuration file.  This is a final process after all steps
        complete and was done to reduce amount of disk space needed to store these packets.  In the subsequent
        sequence pipeline, all molecular packets will populate the flow cell since the wash out has essentially
        already been simulated here.
        :param args: command line arguements
        """
        start = time.time()
        stage_name = "library_prep_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        # Molecule packets from an input folder
        input_directory_path = self.configuration[stage_name]["input"].get("directory_path",None)
        # TODO: allow the specific molecule files to be set explicitly in configuration
        if input_directory_path:
            molecule_packet_file_paths = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}molecule*.txt', recursive=True)
            file_count = len(molecule_packet_file_paths)
            packet_ids = list(range(file_count))
        else:
            molecule_packet_file_paths = []
            file_count = 0
            packet_ids = []

        # Molecule packets from distributions
        from_distribution_data = self.configuration[stage_name]["input"].get("from_distribution_data", None)
        for from_dist in from_distribution_data.values():
            file_count += from_dist['num_packets']

        data_directory = os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME)
        log_directory = os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME)
        directory_structure = GeneralUtils.create_subdirectories(file_count, data_directory)
        self.create_step_log_directories(file_count, stage_name, log_directory)

        self.setup_dispatcher(
                args.dispatcher_mode,
                stage_name,
                os.path.join(self.output_directory_path, stage_name),
                directory_structure
        )
        self.dispatcher.dispatch(
                molecule_packet_file_paths,
                packet_ids,
                from_distribution_data,
        )
        print(f"Finished executing library prep after {time.time() - start:0.2f}s")

    def run_sequence_pipeline(self, args, setup=True):
        """
        This is how run_beers.py calls the sequence pipeline by itself.  Controller attributes are set up.
        All molecule packet file are located - note that these molecule packet files
        are really the outputs of the various library prep processes and as such one can point the input directory for
        the sequence pipeline to the location of the files created by an earlier call to the library pipeline via the
        configuration file.  In that way, the library prep pipeline and the sequence pipeline can be run one after the
        other.  However, here again, the sequence pipeline assumes all molecule packets needed for processing are
        already in place.

        This pipeline stage is more elaborate than that of the library prep pipeline because it handles flowcell
        loading and FASTQ reporting.  Both of these processes are handled by the controller since a knowledge of all
        packets is necessary to perform those functions.

        Here too, the data directory path and the number of packets to be processed are used to fully develop the
        output folder structure, which again may have nested sub-directories in the case of a very large number of
        input molecule packets.  Then each molecule packet is loaded in turn onto the flowcell.  Those molecules to be
        washed out have already been so in the prior pipeline step.  These molecules are located by flowcell coordinates
        and returned as clusters bundled into a cluster packet.  Those cluster packets are then serialized into a gzip
        file located inside the controller data folder in the output directory.  Note that this is really intermediate
        data that will be fed into the sequence pipeline stage proper.

        Finally as with the library prep pipeline, the data and log output directories are identified and the data
        directory and number of packets to be processed are used to fully develop the output folder structure, which
        may have nested sub-directories in the case of a very large number of input cluster packets.  Again, the step
        log directories are created. An auditor is created with a path to an audit file and a list of the number of
        packet processes to expect.  The dispatcher is instantiated and finally, the dispatcher is run with the
        cluster_packet_file_paths provided.  The auditor loops waiting until all sequence pipeline processes are
        complete.  Once that happens, the collected reads are formatted into 1 or 2 FASTQ files for each lane of the
        flowcell used.
        :param args: The command line arguments
        """
        stage_name = "sequence_pipeline"
        if setup:
            self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        intermediate_directory_path = os.path.join(self.output_directory_path, self.controller_name)
        data_directory_path = os.path.join(intermediate_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        cluster_packet_file_paths = []
        molecule_packet_file_paths = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.txt', recursive=True)
        if not molecule_packet_file_paths:
            raise ControllerValidationException(f"No molecule packet files were found in the given input directory, "
                                                f"{input_directory_path}, or any of its subdirectories")
        file_count = len(molecule_packet_file_paths)
        directory_structure = GeneralUtils.create_subdirectories(file_count, data_directory_path)
        self.setup_flowcell()
        for molecule_packet_filename in molecule_packet_file_paths:
            print(f"Loading {molecule_packet_filename} to flowcell")
            molecule_packet = MoleculePacket.deserialize(molecule_packet_filename)
            cluster_packet = self.flowcell.load_flowcell(molecule_packet)
            cluster_packet_filename = f"cluster_packet_start_pkt{cluster_packet.cluster_packet_id}.gzip"
            subdirectory_list = \
                GeneralUtils.get_output_subdirectories(cluster_packet.cluster_packet_id, directory_structure)
            data_subdirectory_path = os.path.join(data_directory_path, *subdirectory_list)
            cluster_packet_file_path = os.path.join(data_subdirectory_path, cluster_packet_filename)
            cluster_packet.serialize(cluster_packet_file_path)
            cluster_packet_file_paths.append(cluster_packet_file_path)
        packet_file_count = len(cluster_packet_file_paths)
        data_directory_path = os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME)
        log_directory_path = os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME)
        directory_structure = GeneralUtils.create_subdirectories(packet_file_count, data_directory_path)
        self.create_step_log_directories(file_count, stage_name, log_directory_path)
        auditor = Auditor(packet_file_count, os.path.join(self.output_directory_path, stage_name))
        self.setup_dispatcher(args.dispatcher_mode,
                              stage_name,
                              os.path.join(self.output_directory_path, stage_name),
                              directory_structure)
        self.dispatcher.dispatch(cluster_packet_file_paths)

        while not auditor.is_processing_complete():
            time.sleep(1)

        fastq_output = self.configuration[stage_name]["output"]["output_fastq"]
        sam_output = self.configuration[stage_name]["output"]["output_sam"] or not fastq_output

        sample_barcode_map = {sample: (cfg['barcodes']['i5'] + "+" + cfg['barcodes']['i7'])
                                    for sample, cfg in self.configuration['samples'].items()}

        if fastq_output:
            print("Generating FastQs")
            for lane in self.flowcell.lanes_to_use:
                fastq_file = os.path.join(
                    self.output_directory_path,
                    self.controller_name,
                    CONSTANTS.DATA_DIRECTORY_NAME,
                )
                fast_q = FastQ(
                        lane,
                        os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME),
                        fastq_file,
                        sample_barcodes = sample_barcode_map,
                )
                fast_q.generate_report()

        if sam_output:
            print("Generating SAMs")
            #TODO: allow BAM creation?
            for lane in self.flowcell.lanes_to_use:
                reference_genome = read_fasta(self.configuration['resources']['reference_genome_fasta'])

                sam_file = os.path.join(self.output_directory_path, self.controller_name, CONSTANTS.DATA_DIRECTORY_NAME)
                sam = SAM(
                    lane,
                    os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME),
                    sam_file,
                    sample_barcodes = sample_barcode_map)
                sam.generate_report(reference_genome, BAM=False)

    def run_prep_and_sequence_pipeline(self, args):
        """
        Run first Library Prep pipeline and then the Sequence Pipeline using the output
        from the library prep as input for sequence.
        :param args: The command line arguments
        """
        self.run_library_prep_pipeline(args)
        # Chain library prep output into sequence input
        self.configuration['sequence_pipeline']['input']['directory_path'] = os.path.join(self.output_directory_path, "library_prep_pipeline", "data")
        self.run_sequence_pipeline(args, setup=False)

    def perform_setup(self, args, stage_names):
        """
        This helper method sets up a number of attributes and behaviors in the controller.  Stacktraces are suppressed
        and only user friendly errors are shown when the debugger is off (just a command line arg right now).  The full
        configuration file data and run id are salted away and the random seed is set.  The initial output folder
        structure (excluding the subdirectory structure needed to accommodate large numbers of files) is created.  The
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
        Validator.validate(stage_names, self.configuration)
        self.plant_seed()
        self.create_output_folder_structure(stage_names, overwrite=args.force_overwrite)
        self.create_controller_log()

    def setup_dispatcher(self, dispatcher_mode, stage_name, output_directory_path, nested_depth):
        """
        The dispatcher is now instantiated and prepared.  If no dispatcher mode (lsf, serial, multicore) is set, it is
        read from the configuration file.  If the mode is not present there, an exception is raised.  The dispatcher
        is then instantiated and and attached to the controller as an attribute.
        :param dispatcher_mode: The type of dispatch to perform (serial, lsf, multicore)
        :param stage_name: The name of the pipeline stage (library_prep_pipeline or sequence_pipeline)
        :param output_directory_path: The top level ouptut directory path
        :param nested_depth: nesting data the pipelines need to write data and log files to their proper locations.
        """
        if not dispatcher_mode:
            if not self.controller_configuration.get('dispatcher_mode', None):
                raise ControllerValidationException('No dispatcher_mode given either on the command line'
                                                    ' or in the configuration file')
            dispatcher_mode = self.controller_configuration['dispatcher_mode']
        if dispatcher_mode not in SUPPORTED_SCHEDULER_MODES:
            raise ControllerValidationException(f'{dispatcher_mode} is not a supported mode.\n'
                                                'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        self.dispatcher = Dispatcher(self.run_id,
                                     dispatcher_mode,
                                     stage_name,
                                     self.configuration,
                                     self.configuration_file_path,
                                     output_directory_path,
                                     nested_depth)

    def setup_flowcell(self):
        """
        Instantiates the flowcell, validates the flowcell parameters and attaches it to the controller.
        """
        self.flowcell = Flowcell(self.run_id, self.configuration, self.configuration[self.controller_name]['flowcell'])
        valid, msg = self.flowcell.validate()
        if not valid:
            raise (ControllerValidationException(msg))

    def retrieve_configuration(self, configuration_file_path):
        """
        Helper method to parse the configuration file given by the path info a dictionary attached to the controller
        object.  For convenience, the portions of the configuration file that contains parametric data specific to the
        controller and the resources are set to separate dictionaries also attached to the controller.
        :param configuration_file_path: The absolute file path of the configuration file
        """
        with open(configuration_file_path, "r+") as configuration_file:
            try:
                self.configuration = json.load(configuration_file)
            except json.decoder.JSONDecodeError:
                print(f"ERROR: JSON Error reading in {configuration_file_path}:")
                raise
        self.configuration_file_path = configuration_file_path
        self.controller_configuration = self.configuration[self.controller_name]
        self.resources = self.configuration[self.resources_name]
        self.resources['resources_folder'] = os.path.join(CONSTANTS.ROOT_DIR, "resources")

    def set_run_id(self, run_id):
        """
        Helper method to add the run id to the controller for later use.  The run id, if any, on the command line
        takes precedence.  If no run id is present on the command line, the controller configuration data is searched
        for it.  If no run id is found in either place, an error is raised.
        :param run_id: The run id found on the command line, if any.
        """
        if not run_id:
            if not self.controller_configuration['run_id']:
                raise ControllerValidationException('No run id given either on the command line'
                                                    ' or in the configuration file')
            self.run_id = self.controller_configuration['run_id']
        else:
            self.run_id = run_id

    def plant_seed(self):
        """
        Helper method to obtain the seed from the controller configuration data or to generate a seed if
        no seed is provided by the user and then set it.  The seed will be added to the controller log file created for
        this run so that the user may re-create the run exactly at a later date, assuming all else remains the same.
        """
        self.seed = self.controller_configuration.get('seed', GeneralUtils.generate_seed())
        np.random.seed(self.seed)

    def create_output_folder_structure(self, stage_names, overwrite=False):
        """
        Use the provided stage names, the run id and the top level output directory path from the configuration data
        to create the top level directory of a preliminary directory structure.  The attempt fails if either the top
        level output directory exists but either is not a directory or is a non-empty directory or if the user has
        insufficient permissions to create the directory. Created in the level directly below the top level output
        directory, are folders named after the stage names provided (i.e., controller, library_prep_pipeline,
        sequence_pipeline) and beneath each of these are data and log folders.  Additional subdirectories are created
        later to organize the numerous files expected and avoid overloading any one directory.

        Example of top level output folder structure:

        lib_prep_results_run101
            controller
                data
                logs
            library_prep_pipline
                data
                logs

        :param stage_names: names of folders directly below the top level output directory (e.g., controller,
        library_prep)
        """
        self.output_directory_path = f"{self.controller_configuration['output_directory_path']}_run{self.run_id}"
        if not os.path.exists(self.output_directory_path):
            try:
                os.makedirs(self.output_directory_path, mode=0o0755, exist_ok=True)
            except PermissionError:
                raise ControllerValidationException(f"Insufficient permissions to create {self.output_directory_path}")
        else:
            if not os.path.isdir(self.output_directory_path):
                raise ControllerValidationException(f"The output directory path, {self.output_directory_path},"
                                                    f" is not a directory.")
            if os.listdir(self.output_directory_path):
                if overwrite:
                    print(f"Clearing output directory {self.output_directory_path}")
                    shutil.rmtree(self.output_directory_path)
                else:
                    raise ControllerValidationException(f"The output directory path, {self.output_directory_path},"
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
        input_directory_path = self.configuration["expression_pipeline"]["input"]["directory_path"]
        self.input_samples = []
        for input_sample in self.configuration['expression_pipeline']["input"]["data"].values():
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

    def create_step_log_directories(self, file_count, stage_name, log_directory_path):
        """
        The number and type of steps in either the library_prep pipeline or the sequence pipeline will depend upon the
        user.  That information is provided in the configuration data.  Since there will be as many step logs as there
        are packets to process, they too must be organized into subdirectories.  We also need the same subdirectory
        organization for the standard out, standard error and pipeline log file.  So that is created here along with
        the subdirectory structure for each step discovered in the configuration data for the given stage.
        :param file_count: number of files to house in the subdirectory structure (i.e., number of pipeline processes
        to run)
        :param stage_name: the stage name of the pipeline process (e.g., library_prep_pipeline, sequence_pipeline)
        :param log_directory_path: The top level log directory to house this sub-structure.
        """
        for subdirectory_name in [LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                  CONSTANTS.STDERR_SUBDIRECTORY_NAME,
                                  CONSTANTS.STDOUT_SUBDIRECTORY_NAME]:
            GeneralUtils.create_subdirectories(file_count, os.path.join(log_directory_path, subdirectory_name))
        for step in self.configuration[stage_name]['steps']:
            step_name = step['step_name'].split(".")[1]
            GeneralUtils.create_subdirectories(file_count, os.path.join(log_directory_path, step_name))


class ControllerValidationException(Exception):
    pass
