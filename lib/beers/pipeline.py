#!/usr/bin/python

import json
import importlib
import numpy as np
from molecule import Molecule

class Pipeline:

    def __init__(self):
        with open("../../config/config.json", "r+") as config_file:
            config_json = json.load(config_file)
            self.steps = []
            for step in config_json["pipeline"]:
                step_name = step["step_name"]
                log_file = step["log_file"]
                parameters = step["parameters"]
                module = importlib.import_module(step_name)
                step_class = getattr(module, step_name)
                self.steps.append(step_class(log_file, parameters))
            self.sample_file = config_json["sample_file"]
            self.sample = self.populate_molecules()
            self.seed = config_json["base_seed"]
            self.log_filename = config_json["log_filename"]
            self.log_sample()
            np.random.seed(self.seed)

    def validate(self, **kwargs):
        if not all([step.validate() for step in self.steps]):
            raise BeersValidationException("Validation error: see stderr for details.")

    def execute(self):
        sample = self.sample
        for step in self.steps:
            sample = step.execute(sample)

    def log_sample(self):
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in self.sample:
                log_file.write(molecule.log_entry())

    def populate_molecules(self):
        molecules = []
        i = 0
        with open(self.sample_file, 'r') as sample_file:
            sequences = sample_file.readlines()
            for text in sequences:
                sequence = text.rstrip()
                cigar = str(len(sequence)) + 'M'
                molecule = Molecule(i, sequence, 1, cigar)
                molecules.append(molecule)
                i += 1
        return molecules

    @staticmethod
    def main():
        pipeline = Pipeline()
        pipeline.validate()
        pipeline.execute()


class BeersValidationException(Exception):
    pass


if __name__ == "__main__":
    Pipeline.main()
