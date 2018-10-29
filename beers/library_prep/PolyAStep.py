from molecule import Molecule
import numpy as np
import sys
from timeit import default_timer as timer


class PolyAStep:

    def __init__(self, log_filename, parameters):
        self.log_filename = log_filename
        self.min_retention_prob = parameters.get("min_retention_prob")
        self.max_retention_prob = parameters.get("max_retention_prob")
        self.length_retention_prob = parameters.get("length_retention_prob")
        self.min_breakage_prob = parameters.get("min_breakage_prob")
        self.max_breakage_prob = parameters.get("max_breakage_prob")
        self.breakpoint_prob = parameters.get("breakpoint_prob")
        print("Poly A selection step instantiated")

    def execute(self, sample):
        """
        Remove nearly all molecules that don't have a poly A tail.  Most molecules with a poly A tail are retained.
        Some bleed over modeled to occur in both directions.  Additionally, retained molecules have a chance of being
        truncated on the 5' end.
        :param sample: rna molecules subject to the poly A selection step
        :return: rna molecules, possibly modified, that remain following the poly A selection step.
        """
        print("Poly A selection step starting")
        retained_sample = []
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in sample:

                # Weighted distribution based on tail length
                tail_length = molecule.poly_a_tail_length()
                retention_odds = min(self.min_retention_prob + self.length_retention_prob * tail_length,
                                     self.max_retention_prob)
                retained = np.random.choice([1, 0], 1, p=[retention_odds, 1 - retention_odds])[0]
                note = ''
                if retained:
                    retained_sample.append(molecule)
                    note += 'retained'
                    note = self.apply_three_prime_bias(molecule, tail_length, note)
                else:
                    note += 'removed'
                log_file.write(molecule.log_entry(note))
        print("Poly A selection step complete")
        return retained_sample

    def apply_three_prime_bias(self, molecule, tail_length, note):
        """
        Model for polyA selection 3' bias.  Assuming that polyA tail plays no role in truncation of 5' end.
        :param molecule: molecule to evaluate for truncation
        :param tail_length: length of polyA tail, which may be 0
        :param note: comment added to log
        :return: note with an addendum if a truncation occurs
        """

        sequence_minus_tail_length = len(molecule.sequence) - tail_length
        breakage_likelihood = min([self.min_breakage_prob * sequence_minus_tail_length, self.max_breakage_prob])
        breakage = np.random.choice([1, 0], 1, p=[breakage_likelihood, 1 - breakage_likelihood])
        if breakage:

            # Geometric distribution of breakpoints - although geometric dist 1 indexed, we don't subtract
            # the 1 because we want to be sure of leaving at least one non-polyA nt.
            breakpoint_from_3p = min(np.random.geometric(self.breakpoint_prob), sequence_minus_tail_length)
            breakpoint_from_5p = sequence_minus_tail_length - 1 - breakpoint_from_3p
            molecule.truncate(breakpoint_from_5p)
            note += ' broken'
        return note

    def validate(self):
        print(f"Poly A step validating parameters")
        if self.min_retention_prob < 0 or self.min_retention_prob > 1:
            print("The minimum retention probability parameter must be between 0 and 1", file=sys.stderr)
            return False
        if self.max_retention_prob < 0 or self.max_retention_prob > 1:
            print("The maximum retention probability parameter must be between 0 and 1", file=sys.stderr)
            return False
        if self.max_retention_prob < self.min_retention_prob:
            print("The maximum retention probability parameter must be >= to the minimum retention probability.",
                  file=sys.stderr)
            return False
        return True


if __name__ == "__main__":
    np.random.seed(100)
    molecules = []
    i = 0
    with open("../../data/molecules.txt", 'r') as sample_file:
        sequences = sample_file.readlines()
        for text in sequences:
            sequence = text.rstrip()
            cigar = str(len(sequence)) + 'M'
            rna_molecule = Molecule(i, sequence, 1, cigar)
            molecules.append(rna_molecule)
            i += 1
    input_data_log_file = "../../data/polyA_step_input_data.log"
    with open(input_data_log_file, "w+") as input_data_log:
        input_data_log.write(Molecule.header)
        for rna_molecule in molecules:
            input_data_log.write(rna_molecule.log_entry())
    output_data_log_file = "../../data/polyA_step_output_data.log"
    input_parameters = {
        "min_retention_prob": 0.01,
        "max_retention_prob": 0.99,
        "length_retention_prob": 0.04,
        "min_breakage_prob": 0.0005,
        "max_breakage_prob": 0.98,
        "breakpoint_prob": 0.0001
      }
    step = PolyAStep(output_data_log_file, input_parameters)
    start = timer()
    step.execute(molecules)
    end = timer()
    print(f"PolyA Selection Step: {end - start} for {len(molecules)} molecules.")
