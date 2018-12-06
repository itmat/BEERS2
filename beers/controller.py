import json
import time
import numpy as np
import termcolor
import os
import sys
import resource
import traceback
import pickle
from datetime import datetime
from beers.utilities.general_utils import GeneralUtils
from beers.expression.expression_pipeline import ExpressionPipeline
from beers.sequence.sequence_pipeline import SequencePipeline
from beers.sample import Sample
from beers.utilities.adapter_generator import AdapterGenerator
from beers.flowcell import Flowcell
from beers.molecule_packet import MoleculePacket
from beers.dispatcher import Dispatcher
from beers.fast_q import FastQ
import glob


class Controller:

    def __init__(self):
        self.controller_name = 'controller'

    def run_expression_pipeline(self, args):
        stage_name = "expression_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        self.assemble_input_samples()
        ExpressionPipeline.main(self.configuration['expression_pipeline'],
                                os.path.join(self.output_directory_path, stage_name),
                                self.input_samples)

    def run_library_prep_pipeline(self, args):
        stage_name = "library_prep_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        self.setup_dispatcher(args.dispatcher_mode, stage_name, input_directory_path, os.path.join(self.output_directory_path, stage_name))
        # Packet files will be in nested directories (at least one level deep)
        molecule_packet_filenames = glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.gzip', recursive=True)
        self.create_subdirectories(len(molecule_packet_filenames), os.path.join(self.output_directory_path, stage_name))
        self.dispatcher.dispatch(molecule_packet_filenames)

    def run_sequence_pipeline(self, args):
        stage_name = "sequence_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        intermediate_directory_path = os.path.join(self.output_directory_path, self.controller_name, "data")
        self.setup_dispatcher(args.dispatcher_mode, stage_name, intermediate_directory_path, os.path.join(self.output_directory_path, stage_name))
        cluster_packet_file_paths = []
        for molecule_packet_filename in glob.glob(f'{input_directory_path}{os.sep}**{os.sep}*.gzip', recursive=True):
            molecule_packet = MoleculePacket.get_serialized_molecule_packet(input_directory_path, molecule_packet_filename)
            cluster_packet = self.setup_flowcell(molecule_packet)
            cluster_packet_filename = f"cluster_packet_start_pk{cluster_packet.cluster_packet_id}.gzip"
            log_subdirectory_path, data_subdirectory_path = \
                Controller.create_subdirectories(cluster_packet.cluster_packet_id, intermediate_directory_path)
            cluster_packet_file_path = os.path.join(data_subdirectory_path, cluster_packet_filename)
            cluster_packet.serialize(cluster_packet_file_path)
            cluster_packet_file_paths.append(cluster_packet_file_path)
        # Molecule packet no longer needed - trying to save RAM
        molecule_packet = None
        self.create_subdirectories(len(cluster_packet_file_paths), os.path.join(self.output_directory_path, stage_name))
        self.dispatcher.dispatch(cluster_packet_file_paths)
        for lane in self.flowcell.lanes_to_use:
            fast_q = FastQ(lane,
                           os.path.join(self.output_directory_path, stage_name, "data"),
                           os.path.join(self.output_directory_path, self.controller_name, "data"))
            fast_q.generate_report()

    def run_prep_and_sequence_pipeline(self, args):
        pass

    def perform_setup(self, args, stage_names):
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

    def setup_dispatcher(self, dispatcher_mode, stage_name, input_directory_path, output_directory_path):
        if not dispatcher_mode:
            if not self.controller_configuration.get('dispatcher_mode', None):
                raise ControllerValidationException('No dispatcher_mode given either on the command line'
                                                    ' or in the configuration file')
            dispatcher_mode = self.controller_configuration['dispatcher_mode']
        self.dispatcher = Dispatcher(dispatcher_mode, stage_name, self.configuration, input_directory_path, output_directory_path)

    def setup_flowcell(self, molecule_packet):
        self.flowcell = Flowcell(self.run_id, self.configuration, self.configuration[self.controller_name]['flowcell'])
        valid, msg = self.flowcell.validate()
        if not valid:
            raise (ControllerValidationException(msg))
        return self.flowcell.load_flowcell(molecule_packet)

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
            os.makedirs(os.path.join(self.output_directory_path, stage_name, 'logs'), mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(self.output_directory_path, stage_name, 'data'), mode=0o0755, exist_ok=True)

    def get_input_from_pickle(self, stage_name):
        # Molecule packets coming from file location named in configuration when not directly from
        # the library pipeline
        input_directory_path = self.configuration[stage_name]["input"]["directory_path"]
        molecule_packet_filename = self.configuration[stage_name]["input"]["molecule_packet_filename"]
        molecule_packet_file_path = os.path.join(input_directory_path, molecule_packet_filename)
        with open(molecule_packet_file_path, 'rb') as molecule_packet_file:
            molecule_packet = pickle.load(molecule_packet_file)
        print(f"{stage_name} input loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")
        return molecule_packet

    @staticmethod
    def create_subdirectories(packet_count, parent_directory_path):
        max_packet_group = packet_count // 200
        log_directory_path = os.path.join(parent_directory_path, "logs")
        data_directory_path = os.path.join(parent_directory_path, 'data')
        for packet_group in range(max_packet_group + 1):
            log_subdirectory_path = os.path.join(log_directory_path, f'pkt_grp{packet_group}')
            data_subdirectory_path = os.path.join(data_directory_path, f'pkt_grp{packet_group}')
            os.makedirs(log_subdirectory_path, mode=0o0755, exist_ok=True)
            os.makedirs(data_subdirectory_path, mode=0o0755, exist_ok=True)
        return log_subdirectory_path, data_subdirectory_path

    def create_controller_log(self):
        log_file_path = os.path.join(self.output_directory_path,'controller','logs','controller.log')
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
                       input_sample.get("gender", None),
                       adapter_generator.get_unique_adapter_labels()))
            Sample.next_sample_id += 1


class ControllerValidationException(Exception):
    pass
