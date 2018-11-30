import json
import importlib
import pickle
import time
import sys
import os
import resource

import numpy as np

from beers.molecule import Molecule


class LibraryPrepPipeline:

    def __init__(self, configuration, output_directory_path, molecule_packet):
        self.molecule_packet = molecule_packet
        self.original_ids = set(str(m.molecule_id) for m in self.molecule_packet.molecules)
        self.print_summary(self.molecule_packet.molecules)
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, "logs")
        data_directory_path = os.path.join(output_directory_path, 'data')
        self.log_file_path = os.path.join(log_directory_path, "library_pipeline.log")
        self.steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_molecule_pkt{self.molecule_packet.molecule_packet_id}.log"
            step_log_file_path = os.path.join(log_directory_path, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package="beers.library_prep")
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters))

        results_filename = f"library_prep_pipeline_result_molecule_pkt{self.molecule_packet.molecule_packet_id}.pickle"
        self.results_file_path = os.path.join(data_directory_path, results_filename)

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
            print(f"Finished pipeline in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a pickle file for inspection
        with open(self.results_file_path, "wb") as results_file:
            pickle.dump(molecule_packet, results_file)
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
    def main(configuration, output_directory_path, molecule_packet):
        library_prep_pipeline = LibraryPrepPipeline(configuration, output_directory_path, molecule_packet)
        library_prep_pipeline.validate()
        library_prep_pipeline.execute()

class BeersValidationException(Exception):
    pass
