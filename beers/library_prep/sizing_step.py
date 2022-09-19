import sys
from timeit import default_timer as timer
from scipy.stats import norm
import numpy as np
#import pylab as pl
import pickle
from beers.utilities.library_prep_utils import Utils

from beers_utils.molecule import Molecule



class SizingStep:
    """
    This step simulates filtering molecules by size.

    Molecules are selected with probability that ranges piecewise linearly from 0
    at 'min_length' up to 1 between 'select_all_start_length' and 'select_all_end_length'
    and then down to 0 again at 'max_length'.

    Probability of retention by length
                     __________________
       1|        ___/                  \___
    p   |   ____/                          \____
       0|__/                                    \___
        ---------------------------------------------
           |          |               |          |
           min_length |               |          max_length
                      |   select_all  |
                    start            end
    If select_all_start_length or select_all_end_length are not provided
    then they are set to min_lenght and max_length respectively and the probability
    jumps from 0 to 1 immediately.
    """

    name = "Sizing Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        """
        Initializes the step with a file path to the step log and a dictionary of parameters.
        :param step_log_file_path: location of step logfile
        :param parameters: dictionary of parameters, which must include at least min_length and max_length
        :param global_config: dictionary of step-independent configuration settings
        """
        self.log_filename = step_log_file_path
        self.min_length = parameters.get("min_length")
        self.max_length = parameters.get("max_length")
        self.select_all_start_length = parameters.get("select_all_start_length", self.min_length)
        self.select_all_end_length = parameters.get("select_all_end_length", self.max_length)
        self.global_config = global_config
        print("Sizing step instantiated")

    def execute(self, molecule_packet, rng):
        """
        Remove those molecules that are outside the filter range.  Parameters dictate whether cutoffs are sharp or
        dictated by a distribution function
        :param molecule_packet: rna molecules subject to sizing
        :return: rna molecules retained following filtration
        """
        print("Sizing step starting")
        retained_molecules = []
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                seq_length = len(molecule.sequence)
                if (seq_length < self.min_length or seq_length > self.max_length):
                    retained = False
                elif (seq_length >= self.select_all_start_length) and (seq_length <= self.select_all_end_length):
                    retained = True
                else:
                    if seq_length < self.select_all_start_length:
                        retention_prob = (seq_length - self.min_length) / (self.select_all_start_length - self.min_length)
                    else:
                        retention_prob = (self.max_length - seq_length) / (self.max_length - self.select_all_end_length)
                    retained = (rng.random() < retention_prob)
                note = ''
                if retained:
                    retained_molecules.append(molecule)
                    note += 'retained'
                else:
                    note += 'removed'
                log_file.write(molecule.log_entry(note))
        print("Sizing step complete")
        molecule_packet.molecules = retained_molecules
        return molecule_packet

    @staticmethod
    def validate(parameters, global_config):

        errors = []
        for param in ['min_length', 'max_length']:
            if param not in parameters:
                errors.append(f"Must specify '{param}'")
            elif not (0 <= parameters[param]):
                errors.append(f"{param} must be positive")
        if len(errors) == 0:
            if parameters['min_length'] >= parameters['max_length']:
                errors.append("min_length must be less than max_length")

            if 'select_all_start_length' in parameters or 'select_all_end_length' in parameters:
                min_length = parameters.get("min_length")
                max_length = parameters.get("max_length")
                select_all_start_length = parameters.get("select_all_start_length", min_length)
                select_all_end_length = parameters.get("select_all_end_length", max_length)
                if not (min_length <= select_all_start_length <= select_all_end_length <= max_length):
                    errors.append("SizingStep needs min_length <= select_all_start_length <= select_all_end_length <= max_length")

        return errors
