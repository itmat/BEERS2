import json
import importlib
import pickle
import time
import sys
import os
import resource
import argparse
import json

import numpy as np

from beers.molecule import Molecule
from beers.molecule_packet import MoleculePacket
from beers.utilities.general_utils import GeneralUtils


class LibraryPrepPipeline:

    stage_name = "library_prep_pipeline"
    package = "beers.library_prep"

    def __init__(self, configuration, output_directory_path, molecule_packet):
        self.molecule_packet = molecule_packet
        self.original_ids = set(str(m.molecule_id) for m in self.molecule_packet.molecules)
        self.print_summary(self.molecule_packet.molecules)
        log_subdirectory_path, data_subdirectory_path = \
            GeneralUtils.create_output_subdirectories(self.molecule_packet.molecule_packet_id, output_directory_path)
        self.log_file_path = os.path.join(log_subdirectory_path,
                                          f"{LibraryPrepPipeline.stage_name}_"
                                          f"molecule_pkt{self.molecule_packet.molecule_packet_id}.log")
        self.steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_molecule_pkt{self.molecule_packet.molecule_packet_id}.log"
            step_log_file_path = os.path.join(log_subdirectory_path, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters))
        results_filename = f"{LibraryPrepPipeline.stage_name}_" \
                           f"result_molecule_pkt{self.molecule_packet.molecule_packet_id}.gzip"
        self.results_file_path = os.path.join(data_subdirectory_path, results_filename)

    def validate(self, **kwargs):
        if not all([molecule.validate() for molecule in self.molecule_packet.molecules]):
            raise BeersValidationException("Validation error in molecule packet: see stderr for details.")
        if not all([step.validate() for step in self.steps]):
            raise BeersValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        with open(self.log_file_path, 'w') as log_file:
            self.log_sample(log_file)
            pipeline_start = time.time()
            molecule_packet = self.molecule_packet
            for step in self.steps:
                step_name = step.__class__.name if hasattr(step.__class__, 'name') else step.__class__.__name__
                step_start = time.time()
                molecule_packet = step.execute(molecule_packet)
                elapsed_time = time.time() - step_start

                self.print_summary(molecule_packet.molecules, elapsed_time)

                random_state = np.random.get_state()
                log_file.write(f"# random state following {step_name} is {random_state}\n")
                print(f"{step_name} complete - process RAM currently at"
                      f" {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")

            pipeline_elapsed_time = time.time() - pipeline_start
            print(f"Finished {LibraryPrepPipeline.stage_name} in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a gzip file for inspection
        molecule_packet.serialize(self.results_file_path)
        print(f"Output final sample to {self.results_file_path}")

    def log_sample(self, log_file):
        log_file.write(Molecule.header)
        for molecule in self.molecule_packet.molecules:
            log_file.write(molecule.log_entry())

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
    def main(configuration, input_directory_path, output_directory_path, molecule_packet_filename):
        configuration = json.loads(configuration)
        molecule_packet = MoleculePacket.get_serialized_molecule_packet(input_directory_path, molecule_packet_filename)
        library_prep_pipeline = LibraryPrepPipeline(configuration, output_directory_path, molecule_packet)
        library_prep_pipeline.validate()
        library_prep_pipeline.execute()

class BeersValidationException(Exception):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Library Prep Pipeline')
    parser.add_argument('-c', '--config', required=True, help='Configuration')
    parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
    parser.add_argument('-p', '--molecule_packet_filename', required=True, help="Serialized Molecule Packet Filename.")
    args = parser.parse_args()
    LibraryPrepPipeline.main(args.config, args.output_directory, args.molecule_packet_filename)

