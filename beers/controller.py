import json
import time
import numpy as np
import termcolor
import os
import sys
import traceback
from datetime import datetime
from beers.constants import CONSTANTS
from beers.utilities.general_utils import GeneralUtils
from beers.expression.expression_pipeline import ExpressionPipeline
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sample import Sample
from beers.utilities.adapter_generator import AdapterGenerator
from beers.flowcell import Flowcell
from beers.molecule_packet import MoleculePacket
from beers.dispatcher import Dispatcher
from beers.fast_q import FastQ
from beers.auditor import Auditor
import glob


class Controller:
    """
    The object essentially controls the flow of the pipeline.  The run_beers.py command instantiates a controller and
    calls one of the controller methods depending upon the pipeline-stage requested.  The other methods in the class
    are helper methods.
    """

    def __init__(self):
        """
        These attributes have no value at the time of instantiation.  They are populated by the helper methods.
        """
        self.controller_name = 'controller'
        self.controller_log_filename = 'controller.log'
        # The following attributes are defined following instantiation
        self.debug = False
        self.run_id = None
        self.dispatcher = None
        self.configuration = None
        self.controller_configuration = None
        self.flowcell = None
        self.seed = None
        self.output_directory_path = None
        self.input_samples = []

    def run_expression_pipeline(self, args):
        """
        This is how run_beers.py calls the expression pipeline by itself.  Little is in place here because the
        expression pipeline is in a state of 'undress'.
        :param args: command line arguments
        """
        stage_name = "expression_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        self.assemble_input_samples()
        ExpressionPipeline.main(self.configuration['expression_pipeline'],
                                os.path.join(self.output_directory_path, stage_name),
                                self.input_samples)

    def run_library_prep_pipeline(self, args):
        """
        This is how run_beers.py calls the library prep pipeline by itself.  This pipeline is complete and functional.
        Controller attributes are set up.  All molecule packet files are located.  The data and log output directories
        are identified and the data directory and number of packets to be processed are used to fully develop the
        output folder structure, which may have nested sub-directories in the case of a very large number of input
        molecule packets.  The step log directories are created.  The dispatcher is instantiated and finally, the
        dispatcher is run with the molecule_packet_file_paths provided.  Note that this pipeline stage assumes all
        molecule packets are immediately available.  That will likely not be the case when the expression and library
        prep pipelines are run together.  That will be a different run_beers.py call.
        :param args: command line arguements
        """
        stage_name = "library_prep_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        molecule_packet_file_paths = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.gzip', recursive=True)
        file_count = len(molecule_packet_file_paths)
        data_directory = os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME)
        log_directory = os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME)
        directory_structure = GeneralUtils.create_subdirectories(file_count, data_directory)
        self.create_step_log_directories(file_count, stage_name, log_directory)
        self.setup_dispatcher(args.dispatcher_mode,
                              stage_name,
                              input_directory_path,
                              os.path.join(self.output_directory_path, stage_name),
                              directory_structure)
        self.dispatcher.dispatch(molecule_packet_file_paths)

    def run_sequence_pipeline(self, args):
        """
        This is how run_beers.py calls the sequence pipeline by itself.  This pipeline is complete and functional.
        Controller attributes are set up.  All molecule packet file are located - note that these molecule packet files
        are really the outputs of the various library prep processes as such one can point the input directory for the
        sequence pipeline to the location of the files created by an earlier call to the library pipeline via the
        configuration file.  In that way, the library prep pipeline and the sequence pipeline can be run one after the
        other.  However, here again, the sequence pipeline assumes all molecule packets needed for processing are
        already in place.  That will not necessarily be the case when the library_prep_and_sequence pipeline is called
        as the molecule packets and subsequent cluster packets may be handled individually.

        This pipeline stage is more elaborate than that of the library prep pipeline because it handles flowcell
        loading and FASTQ reporting.  Both of these processes are handled by the controller since a knowledge of all
        packets is necessary to perform those functions.

        Here too, the data directory path and the number of packets to be processed are used to fully develop the
        output folder structure, which again may have nested sub-directories in the case of a very large number of
        input molecule packets.  Then each molecule packet is loaded in turn onto the flowcell.  Those molecules that
        are retained by the flowcell are located by flowcell coordinates and returned as clusters bundled into a
        cluster packet.  Those cluster packets are then serialized into a gzip file located inside the controller
        data folder in the output directory.  Note that this is really intermediate data that will be fed into the
        sequence pipeline stage proper.

        Finally as with the library prep pipeline, the data and log output directories are identified and the data
        directory and number of packets to be processed are used to fully develop the output folder structure, which
        may have nested sub-directories in the case of a very large number of input cluster packets.  Again, the step
        log directories are created. An auditor is created with a path to an audit file and a list of the number of
        packet processes to expect.  The dispatcher is instantiated and finally, the dispatcher is run with the
        cluster_packet_file_paths provided.  The auditor loops waiting under all sequence pipeline processes are
        complete.  Once that happens the collected reads are formatted into 1 or 2 FASTQ files for each lane of the
        flowcell used.
        :param args: The command line arguments
        """
        stage_name = "sequence_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        intermediate_directory_path = os.path.join(self.output_directory_path, self.controller_name)
        data_directory_path = os.path.join(intermediate_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        cluster_packet_file_paths = []
        molecule_packet_file_paths = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.gzip', recursive=True)
        if not molecule_packet_file_paths:
            raise ControllerValidationException(f"No molecule packet files were found in the given input directory, "
                                                f"{input_directory_path}, or any of its subdirectories")
        file_count = len(molecule_packet_file_paths)
        directory_structure = GeneralUtils.create_subdirectories(file_count, data_directory_path)
        self.setup_flowcell()
        for molecule_packet_filename in molecule_packet_file_paths:
            molecule_packet = MoleculePacket.get_serialized_molecule_packet(input_directory_path,
                                                                            molecule_packet_filename)
            cluster_packet = self.flowcell.load_flowcell(molecule_packet)
            cluster_packet_filename = f"cluster_packet_start_pkt{cluster_packet.cluster_packet_id}.gzip"
            subdirectory_list = \
                GeneralUtils.get_output_subdirectories(cluster_packet.cluster_packet_id, directory_structure)
            data_subdirectory_path = os.path.join(data_directory_path, *subdirectory_list)
            cluster_packet_file_path = os.path.join(data_subdirectory_path, cluster_packet_filename)
            cluster_packet.serialize(cluster_packet_file_path)
            cluster_packet_file_paths.append(cluster_packet_file_path)
        # Molecule packet no longer needed - trying to save RAM
        molecule_packet = None
        packet_file_count = len(cluster_packet_file_paths)
        data_directory_path = os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME)
        log_directory_path = os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME)
        directory_structure = GeneralUtils.create_subdirectories(packet_file_count, data_directory_path)
        self.create_step_log_directories(file_count, stage_name, log_directory_path)
        auditor = Auditor(packet_file_count, os.path.join(self.output_directory_path, stage_name))
        self.setup_dispatcher(args.dispatcher_mode,
                              stage_name,
                              intermediate_directory_path,
                              os.path.join(self.output_directory_path, stage_name),
                              directory_structure)
        self.dispatcher.dispatch(cluster_packet_file_paths)
        while not auditor.is_processing_complete():
            time.sleep(1)
        for lane in self.flowcell.lanes_to_use:
            fast_q = FastQ(lane,
                           os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME),
                           os.path.join(self.output_directory_path,
                                        self.controller_name,
                                        CONSTANTS.DATA_DIRECTORY_NAME))
            fast_q.generate_report()

    def run_prep_and_sequence_pipeline(self, args):
        pass

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
        self.plant_seed()
        self.create_output_folder_structure(stage_names)
        self.create_controller_log()

    def setup_dispatcher(self, dispatcher_mode, stage_name, input_directory_path, output_directory_path, nested_depth):
        """
        The dispatcher is now instantiated and prepared.  If no dispatcher mode (lsf, serial, multicore) is set, it is
        read from the configuration file.  If the mode is not present there, an exception is raised.  The dispatcher
        is then instantiated and and attached to the controller as an attribute.
        :param dispatcher_mode: The type of dispatch to perform (serial, lsf, multicore)
        :param stage_name: The name of the pipeline stage (library_prep_pipeline or sequence_pipeline)
        :param input_directory_path: The top level directory where the input packets are found
        :param output_directory_path: The top level ouptut directory path
        :param nested_depth: nesting data the pipelines need to write data and log files to their proper locations.
        """
        if not dispatcher_mode:
            if not self.controller_configuration.get('dispatcher_mode', None):
                raise ControllerValidationException('No dispatcher_mode given either on the command line'
                                                    ' or in the configuration file')
            dispatcher_mode = self.controller_configuration['dispatcher_mode']
        self.dispatcher = Dispatcher(self.run_id,
                                     dispatcher_mode,
                                     self.seed,
                                     stage_name,
                                     self.configuration,
                                     input_directory_path,
                                     output_directory_path,
                                     nested_depth)

    def setup_flowcell(self):
        """
        Instantiates the flowcell, valudates the flowcell parameters and attaches it to the controller.
        """
        self.flowcell = Flowcell(self.run_id, self.configuration, self.configuration[self.controller_name]['flowcell'])
        valid, msg = self.flowcell.validate()
        if not valid:
            raise (ControllerValidationException(msg))

    def retrieve_configuration(self, configuration_file_path):
        with open(configuration_file_path, "r+") as configuration_file:
            self.configuration = json.load(configuration_file)
        self.controller_configuration = self.configuration[self.controller_name]

    def set_run_id(self, run_id):
        if not run_id:
            if not self.controller_configuration['run_id']:
                raise ControllerValidationException('No run id given either on the command line'
                                                    ' or in the configuration file')
            self.run_id = self.controller_configuration['run_id']
        else:
            self.run_id = run_id

    def plant_seed(self):
        self.seed = self.controller_configuration.get('seed', GeneralUtils.generate_seed())
        np.random.seed(self.seed)

    def create_output_folder_structure(self, stage_names):
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
                raise ControllerValidationException(f"The output directory path, {self.output_directory_path},"
                                                    f" must be empty.")
        for stage_name in stage_names:
            os.makedirs(os.path.join(self.output_directory_path, stage_name, CONSTANTS.LOG_DIRECTORY_NAME),
                        mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(self.output_directory_path, stage_name, CONSTANTS.DATA_DIRECTORY_NAME),
                        mode=0o0755, exist_ok=True)

    def create_controller_log(self):
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

    def assemble_input_samples(self):
        # TODO remove hardcoded file - probably a filename lookup based on prep kit used as provided by user
        adapter_generator = AdapterGenerator("TruSeq_adapter_sequences_with_barcodes.MiSeq_HiSeq2000_HiSeq2500.fa")
        input_directory_path = self.controller_configuration["input"]["directory_path"]
        self.input_samples = []
        for input_sample in self.controller_configuration["input"]["data"]:
            sample_name = os.path.splitext(input_sample["filename"])[0]
            input_sample_file_path = os.path.join(input_directory_path, input_sample["filename"])
            self.input_samples.append(
                Sample(Sample.next_sample_id,
                       sample_name,
                       input_sample_file_path,
                       adapter_generator.get_unique_adapter_sequences(),
                       input_sample.get("gender", None)))
            Sample.next_sample_id += 1

    def create_step_log_directories(self, file_count, stage_name, log_directory_path):
        for subdirectory_name in [LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                  CONSTANTS.STDERR_SUBDIRECTORY_NAME,
                                  CONSTANTS.STDOUT_SUBDIRECTORY_NAME]:
            GeneralUtils.create_subdirectories(file_count, os.path.join(log_directory_path, subdirectory_name))
        for step in self.configuration[stage_name]['steps']:
            step_name = step['step_name'].split(".")[1]
            GeneralUtils.create_subdirectories(file_count, os.path.join(log_directory_path, step_name))


class ControllerValidationException(Exception):
    pass
