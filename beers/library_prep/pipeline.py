#!/usr/bin/python

import json
import importlib
import pickle
import time

import numpy as np

from .molecule import Molecule

class Pipeline:

    def __init__(self, config_filename = "../../config/config.json"):
        with open(config_filename, "r+") as config_file:
            config_json = json.load(config_file)
            self.steps = []
            for step in config_json["library_prep_pipeline"]:
                step_name = step["step_name"]
                log_file = step["log_file"]
                parameters = step["parameters"]
                module = importlib.import_module(f'.{step_name}', package="beers.library_prep")
                step_class = getattr(module, step_name)
                self.steps.append(step_class(log_file, parameters))
            self.sample_file = config_json["sample_file"]
            self.sample = self.populate_molecules()
            self.seed = config_json["seed"]
            self.results_filename = config_json["results_file"]
            self.log_filename = config_json["log_filename"]
            self.log_sample()
            np.random.seed(self.seed)

    def validate(self, **kwargs):
        if not all([molecule.validate() for molecule in self.sample]):
            raise BeersValidationException("Validation error in sample: see stderr for details.")
        if not all([step.validate() for step in self.steps]):
            raise BeersValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        pipeline_start = time.time()
        sample = self.sample
        for step in self.steps:
            step_start = time.time()
            sample = step.execute(sample)
            elapsed_time = time.time() - step_start

            self.print_summary(sample, elapsed_time)

        pipeline_elapsed_time = time.time() - pipeline_start
        print(f"Finished pipeline in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a pickle file for inspection
        with open(self.results_filename, "wb") as results_file:
            pickle.dump(sample, results_file)
        print(f"Output final sample to {self.results_filename}")

    def log_sample(self):
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in self.sample:
                log_file.write(molecule.log_entry())

    def populate_molecules(self):
        molecules = []
        i = 0
        if self.sample_file.endswith(".pickle"):
            print(f"Reading molecule list from pickle file {self.sample_file}")
            with open(self.sample_file, 'rb') as sample_file:
                molecules = list(pickle.load(sample_file))
        else:
            print(f"Reading molecule list from text file {self.sample_file}")
            with open(self.sample_file, 'r') as sample_file:
                sequences = sample_file.readlines()
                for text in sequences:
                    sequence = text.rstrip()
                    cigar = str(len(sequence)) + 'M'
                    molecule = Molecule(i, sequence, 1, cigar)
                    molecules.append(molecule)
                    i += 1

        self.original_ids = set(str(m.molecule_id) for m in molecules)

        self.print_summary(molecules)
        return molecules

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
            print(f" <{size_bin_cutoffs[i]}: {size_counts[i]}")
        print(f">={size_bin_cutoffs[-1]}: {size_counts[-1]}")

    @staticmethod
    def main():
        pipeline = Pipeline()
        pipeline.validate()
        pipeline.execute()


class BeersValidationException(Exception):
    pass


if __name__ == "__main__":
    Pipeline.main()
