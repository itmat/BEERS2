import sys
import pickle
from timeit import default_timer as timer
import numpy as np
from beers_utils.molecule import Molecule

class PolyAStep:
    """
    This step simulates the polyA selection step intended to separate mRNA from all other RNA.  It includes biases but
    idealized behavior is available by not specifying parameter values and the default values provide idealized
    behavior.  Any parameter not specified via the configuration file will default to its idealized setting.
    """

    name = "PolyA Selection Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        """
        Initializes the step with a file path to the step log and a dictionary of parameters.  Idenalized defaults
        substitute for missing parameters.
        :param step_log_file_path: location of step logfile
        :param parameters: dictionary of parameters, all of which are optional.
        :param global_config: dictionary of step-independent configuration settings
        """
        self.log_filename = step_log_file_path
        self.min_polya_tail_length = parameters.get("min_polya_tail_length", 40)
        self.min_retention_prob = parameters.get("min_retention_prob", 0.0)
        self.max_retention_prob = parameters.get("max_retention_prob", 1.0)
        self.length_retention_prob = parameters.get("length_retention_prob", 1.0)
        self.breakpoint_prob_per_base = parameters.get("breakpoint_prob_per_base", 0.0)
        self.global_config = global_config
        print("Poly A selection step instantiated")

    def execute(self, molecule_packet, rng):
        """
        Remove nearly all molecules that don't have a poly A tail.  Most molecules with a poly A tail are retained.
        Some bleed over modeled to occur in both directions.  Additionally, retained molecules have a chance of being
        truncated on the 5' end.
        :param molecule_packet: rna molecules subject to the poly A selection step
        :return: molecule packet containing rna molecules, possibly modified, that remain following the poly A
         selection step.
        """
        print("Poly A selection step starting")
        retained_molecules = []
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:

                # Weighted distribution based on tail length
                tail_length = molecule.poly_a_tail_length()
                tail_length = tail_length if tail_length > self.min_polya_tail_length else 0
                retention_odds = min(self.min_retention_prob + self.length_retention_prob * tail_length,
                                     self.max_retention_prob)
                retained = rng.random() <= retention_odds
                note = ''
                if retained:
                    retained_molecules.append(molecule)
                    note += 'retained'
                    note = self.apply_three_prime_bias(molecule, tail_length, note, rng)
                else:
                    note += 'removed'
                log_file.write(molecule.log_entry(note))
        print("Poly A selection step complete")
        molecule_packet.molecules = retained_molecules
        return molecule_packet

    def apply_three_prime_bias(self, molecule, tail_length, note, rng):
        """
        Model for polyA selection 3' bias.  Assuming that polyA tail plays no role in truncation of 5' end.
        :param molecule: molecule to evaluate for truncation
        :param tail_length: length of polyA tail, which may be 0
        :param note: comment added to log
        :return: note with an addendum if a truncation occurs
        """

        sequence_minus_tail_length = len(molecule.sequence) - tail_length
        # Determine the base at which the molecule would break
        if self.breakpoint_prob_per_base > 0:
            break_point = rng.geometric(self.breakpoint_prob_per_base)
            if break_point < sequence_minus_tail_length:
                # Break occured before the end of the molecule
                # breakpoint is referenced from the end of the pre-polyA tail molecule
                # so we truncate all that occur before it
                breakpoint_from_5p = sequence_minus_tail_length - break_point + 1
                molecule.truncate(breakpoint_from_5p) 
                note += ' broken'
        return note

    @staticmethod
    def validate(parameters, global_config):
        """
        Insures that the parameters provided are valid.  Error messages are sent to stderr.
        :return: True if the step's parameters are all valid and false otherwise.
        """
        errors = []
        min_retention_prob = parameters.get("min_retention_prob", 0.0)
        max_retention_prob = parameters.get("max_retention_prob", 1.0)
        length_retention_prob = parameters.get("length_retention_prob", 1.0)
        breakpoint_prob_per_base = parameters.get("breakpoint_prob_per_base", 0.0)
        min_polya_tail_length = parameters.get("min_polya_tail_length", 40)
        if min_retention_prob < 0 or min_retention_prob > 1:
            errors.append("The minimum retention probability parameter must be between 0 and 1")
        if max_retention_prob < 0 or max_retention_prob > 1:
            errors.append("The maximum retention probability parameter must be between 0 and 1")
        if max_retention_prob < min_retention_prob:
            errors.append("The maximum retention probability parameter must be >= to the minimum retention probability.")
        return errors
