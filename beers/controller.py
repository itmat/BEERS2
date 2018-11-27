import json
import time
from timeit import default_timer as timer
import numpy as np
import os
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
        self.molecule_packet = None
        self.cluster_packet = None

    def run_expression_pipeline(self, args):
        self.perform_setup(args)
        self.assemble_input_samples()
        ExpressionPipeline.main(self.input_samples, self.configuration['expression_pipeline'])

    def run_library_prep_pipeline(self, args):
        self.perform_setup(args)
        LibraryPrepPipeline.main(self.configuration['library_prep_pipeline'])

    def run_sequence_pipeline(self, args):
        self.perform_setup(args)
        if not self.molecule_packet:
            self.get_sequencer_input_from_pickle()
        cluster_packet = self.setup_flowcell()
        SequencePipeline.main(self.configuration['sequence_pipeline'], cluster_packet)

    def perform_setup(self, args):
        self.retrieve_configuration(args.config)
        self.set_run_id(args.run_id)
        self.plant_seed()
        self.create_controller_log()

    def setup_flowcell(self):
        flowcell = Flowcell(self.run_id, self.configuration, self.configuration['controller']['flowcell'])
        valid, msg = flowcell.validate()
        if not valid:
            raise (ControllerValidationException(msg))
        return flowcell.load_flowcell(self.molecule_packet)

    def retrieve_configuration(self, configuration_file_path):
        with open(configuration_file_path, "r+") as configuration_file:
            self.configuration = json.load(configuration_file)

    def set_run_id(self, run_id):
        if not run_id:
            if not self.configuration['controller']['run_id']:
                raise ControllerValidationException('No run id given either on the command line'
                                                    ' or in the configuration file')
            self.run_id = self.configuration['controller']['run_id']
        else:
            self.run_id = run_id

    def plant_seed(self):
        self.seed = self.configuration['controller'].get('seed', GeneralUtils.generate_seed())
        np.random.seed(self.seed)

    def get_sequencer_input_from_pickle(self):
        # Molecule packets coming from file location named in configuration when not directly from
        # the library pipeline
        input_directory_path = self.configuration["sequence_pipeline"]["input"]["directory_path"]
        molecule_packet_filename = self.configuration["sequence_pipeline"]["input"]["molecule_packet_filename"]
        molecule_packet_file_path = os.path.join(input_directory_path, molecule_packet_filename)
        with open(molecule_packet_file_path, 'rb') as molecule_packet_file:
            self.molecule_packet = pickle.load(molecule_packet_file)

    def create_controller_log(self):
        log_file_path = self.configuration['controller']['log_file_path']
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
        input_directory_path = self.configuration['controller']["input"]["directory_path"]
        self.input_samples = []
        for input_sample in self.configuration['controller']["input"]["data"]:
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



