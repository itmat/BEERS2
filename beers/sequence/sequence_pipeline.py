import json
import importlib
import pickle
import time

import numpy as np

from beers.molecule import Molecule

class SequencePipeline:

    def __init__(self, config_filename = "../../config/config.json"):
        with open(config_filename, "r+") as config_file:
            config_json = json.load(config_file)
            self.steps = []
            for step in config_json["sequence_pipeline"]:
                module_name, step_name = step["step_name"].rsplit(".")
                log_file = step["log_file"]
                parameters = step["parameters"]
                module = importlib.import_module(f'.{module_name}', package="beers.library_prep")
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
            raise BeersSequenceValidationException("Validation error in sample: see stderr for details.")
        if not all([step.validate() for step in self.steps]):
            raise BeersSequenceValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        pipeline_start = time.time()
        sample = self.sample
        for step in self.steps:
            step_start = time.time()
            sample = step.execute(sample)
            elapsed_time = time.time() - step_start



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

        return molecules


    @staticmethod
    def main():
        sequence_pipeline = SequencePipeline()
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass


if __name__ == "__main__":
    SequencePipeline.main()
