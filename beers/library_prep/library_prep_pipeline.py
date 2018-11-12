import json
import importlib
import pickle
import time
import sys
import os

import numpy as np

from beers.utilities.general_utils import GeneralUtils
from beers.molecule_packet import MoleculePacket
from beers.molecule import Molecule


class LibraryPrepPipeline:

    def __init__(self, configuration, pipeline_log_file, molecule_packet=None):
            self.steps = []
            library_prep_configuration = configuration['library_prep_pipeline']
            input_directory_path = library_prep_configuration["input"]["directory_path"]
            output_directory_path = library_prep_configuration["output"]["directory_path"]
            self.pipeline_log_file = pipeline_log_file
            for step in library_prep_configuration['steps']:
                module_name, step_name = step["step_name"].rsplit(".")
                step_log_filename = step["log_filename"]
                step_log_file_path = os.path.join(output_directory_path, step_log_filename)
                parameters = step["parameters"]
                module = importlib.import_module(f'.{module_name}', package="beers.library_prep")
                step_class = getattr(module, step_name)
                self.steps.append(step_class(step_log_file_path, parameters))
            self.molecule_packet = molecule_packet
            if not self.molecule_packet:
                molecule_packet_filename = library_prep_configuration["input"]["molecule_packet_filename"]
                self.molecule_packet_file_path = os.path.join(input_directory_path, molecule_packet_filename)
                self.populate_molecule_packet()
            results_filename = library_prep_configuration["output"]["results_filename"]
            self.results_file_path = os.path.join(output_directory_path, results_filename)
            self.log_sample()

    def validate(self, **kwargs):
        if not all([molecule.validate() for molecule in self.molecule_packet.molecules]):
            raise BeersValidationException("Validation error in molecule packet: see stderr for details.")
        if not all([step.validate() for step in self.steps]):
            raise BeersValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        pipeline_start = time.time()
        molecule_packet = self.molecule_packet
        for step in self.steps:
            step_start = time.time()
            molecule_packet = step.execute(molecule_packet)
            elapsed_time = time.time() - step_start

            self.print_summary(molecule_packet.molecules, elapsed_time)

            random_state = np.random.get_state()
            self.pipeline_log_file.write(f"# random state following {step.__class__.__name__} is {random_state}\n")

        pipeline_elapsed_time = time.time() - pipeline_start
        print(f"Finished pipeline in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a pickle file for inspection
        with open(self.results_file_path, "wb") as results_file:
            pickle.dump(molecule_packet, results_file)
        print(f"Output final sample to {self.results_file_path}")

    def log_sample(self):
        self.pipeline_log_file.write(Molecule.header)
        for molecule in self.molecule_packet.molecules:
            self.pipeline_log_file.write(molecule.log_entry())

    def populate_molecule_packet(self):
        print(f"Reading molecule packet from pickle file {self.molecule_packet_file_path}")
        with open(self.molecule_packet_file_path, 'rb') as molecule_packet_file:
            self.molecule_packet = pickle.load(molecule_packet_file)
        self.original_ids = set(str(m.molecule_id) for m in self.molecule_packet.molecules)
        self.print_summary(self.molecule_packet.molecules)

    def print_summary(self, sample, elapsed_time=None):
        '''Output a summary of the sample (number of molecules, time taken, etc.)'''
        if elapsed_time is not None:
            print(f"Step took {elapsed_time:.3} seconds")

        print(f"Sample has {len(sample)} molecules")
        parent_ids = set(str(m.molecule_id).split(".")[0] for m in sample)
        percent_original_represented = len(self.original_ids.intersection(parent_ids))/len(self.original_ids)
        print(f"Percent of the original ids that are still represented: {percent_original_represented:0.2%}")

        size_bin_cutoffs = [100,500,1000]
        size_counts = [0]*(len(size_bin_cutoffs)+1)
        for molecule in sample:
            size_counts[np.searchsorted(size_bin_cutoffs, len(molecule))] += 1
        print(f"Counts of molecules in size ranges:")
        for i in range(len(size_bin_cutoffs)):
            print(f" <={size_bin_cutoffs[i]}: {size_counts[i]}")
        print(f">{size_bin_cutoffs[-1]}: {size_counts[-1]}")

    @staticmethod
    def main(configuration_file_path):
        with open(configuration_file_path, "r+") as configuration_file:
            configuration = json.load(configuration_file)
        seed = configuration.get('seed', GeneralUtils.generate_seed())
        np.random.seed(seed)
        output_directory_path = configuration['library_prep_pipeline']["output"]["directory_path"]
        log_file_path = os.path.join(output_directory_path, configuration['library_prep_pipeline']["output"]["log_filename"])
        with open(log_file_path, "w") as pipeline_log_file:
            pipeline_log_file.write(f"#Seed: {seed} \n")
            library_prep_pipeline = LibraryPrepPipeline(configuration, pipeline_log_file)
            library_prep_pipeline.validate()
            library_prep_pipeline.execute()

class BeersValidationException(Exception):
    pass


if __name__ == "__main__":
    sys.exit(LibraryPrepPipeline.main("../../config/config.json"))
