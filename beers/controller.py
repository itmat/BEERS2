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
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline
from beers.sample import Sample
from beers.utilities.adapter_generator import AdapterGenerator
from beers.flowcell import Flowcell


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

    def run_library_prep_pipeline(self, args, molecule_packet=None):
        stage_name = "library_prep_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        if not molecule_packet:
            molecule_packet = self.get_input_from_pickle(stage_name)
        LibraryPrepPipeline.main(self.configuration[stage_name],
                                 os.path.join(self.output_directory_path, stage_name),
                                 molecule_packet)

    def run_sequence_pipeline(self, args, molecule_packet=None):
        stage_name = "sequence_pipeline"
        self.perform_setup(args, [self.controller_name, stage_name])
        if not molecule_packet:
            molecule_packet = self.get_input_from_pickle(stage_name)
        cluster_packet = self.setup_flowcell(molecule_packet)
        # Molecule packet no longer needed - trying to save RAM
        molecule_packet = None
        SequencePipeline.main(self.configuration[stage_name],
                              os.path.join(self.output_directory_path, stage_name),
                              cluster_packet)

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

    def setup_flowcell(self, molecule_packet):
        flowcell = Flowcell(self.run_id, self.configuration, self.configuration[self.controller_name]['flowcell'])
        valid, msg = flowcell.validate()
        if not valid:
            raise (ControllerValidationException(msg))
        return flowcell.load_flowcell(molecule_packet)

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
