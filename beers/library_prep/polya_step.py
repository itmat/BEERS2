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

    Configuration Example::

        # Chance per base of fragmentation prior to selection
        # Increasing this above 0 induces a 3' bias.
        # A value of 0.001 induces a reasonably high 3' bias
        breakpoint_prob_per_base: 0.0

        # Probability of retention is computed as a value between
        # min_retention_prob and max_retention_prob
        # For every base of the polyA tail beyond min_polya_tail_length
        # the probability increases linearly from min_retention_prob
        # by length_retention_prob, up to a max of max_polya_tail_length.
        max_retention_prob: 1.0
        min_retention_prob: 0.0
        min_polya_tail_length: 40
        length_retention_prob: 0.05
    """

    name = "PolyA Selection Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        self.log_filename = step_log_file_path
        self.min_polya_tail_length = parameters.get("min_polya_tail_length", 40)
        self.min_retention_prob = parameters.get("min_retention_prob", 0.0)
        self.max_retention_prob = parameters.get("max_retention_prob", 1.0)
        self.length_retention_prob = parameters.get("length_retention_prob", 1.0)
        self.breakpoint_prob_per_base = parameters.get("breakpoint_prob_per_base", 0.0)
        self.global_config = global_config
        print("Poly A selection step instantiated")

    def execute(self, molecule_packet, rng):
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

    def apply_three_prime_bias(
            self,
            molecule: Molecule,
            tail_length: int,
            note: str,
            rng: np.random.Generator,
        ) -> str:
        """
        Model for polyA selection 3' bias.  Assuming that polyA tail plays no role in truncation of 5' end.

        Parameters
        ----------
        molecule:
            molecule to evaluate for truncation
        tail_length:
            length of polyA tail, which may be 0
        note:
            comment added to log
        rng:
            random number generator

        Returns
        -------
        str
            note to add to log
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
        if not (0 <= breakpoint_prob_per_base <= 1):
            errors.append("breakpoint_prob_per_base must be between 0 and 1")
        if not (0 <= length_retention_prob <= 1):
            errors.append("length_retention_prob must be between 0 and 1")

        return errors
