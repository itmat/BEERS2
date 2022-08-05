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
        self.idealized = False
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

    def validate(self):
        """
        Insures that the parameters provided are valid.  Error messages are sent to stderr.
        :return: True if the step's parameters are all valid and false otherwise.
        """
        print(f"Sizing step validating parameters")
        if self.idealized:
            if not self.min_length or not self.max_length:
                print("The minimum and maximum cutoff lengths must be specified.", file=sys.stderr)
                return False
            if self.min_length < 0 or self.min_length > self.max_length:
                print("The minimum cutoff length {self.min_length} must be non-zero and less than the"
                      "maximum cutoff length, {self.max_length}.", file=sys.stderr)
                return False
            if not (self.min_length <= self.select_all_start_length <= self.select_all_end_length <= self.max_length):
                print("SizingStep needs min_length <= select_all_start_length <= select_all_end_length <= max_length")
                return False
        else:
            pass
        return True
