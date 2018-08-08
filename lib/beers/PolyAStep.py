from molecule import Molecule
import numpy as np
import sys


class PolyAStep:

    def __init__(self, log_file, parameters):
        self.history_filename = log_file
        self.min_retention_prob = parameters.get("min_retention_prob")
        self.max_retention_prob = parameters.get("max_retention_prob")
        self.length_retention_prob = parameters.get("length_retention_prob")
        self.min_breakage_prob = parameters.get("min_breakage_prob")
        self.max_breakage_prob = parameters.get("max_breakage_prob")
        self.breakpoint_prob = parameters.get("breakpoint_prob")
        print("Poly A selection step instantiated")

    def execute(self, sample):
        print("Poly A selection step starting")
        retained_sample = []
        with open(self.history_filename, "w+") as history_file:
            history_file.write(Molecule.header)
            print(f"sample size in {len(sample)}")
            for molecule in sample:
                tail_length = molecule.poly_a_tail_length()
                retention_odds = min(self.min_retention_prob + self.length_retention_prob * tail_length, self.max_retention_prob)
                retained = np.random.choice([1, 0], 1, p=[retention_odds, 1 - retention_odds])[0]
                note = ''
                if retained:
                    retained_sample.append(molecule)
                    note += 'retained'
                    note = self.three_prime_bias(molecule, tail_length, note)
                else:
                    note += 'removed'
                #print(tail_length, retention_odds, note)
                history_file.write(molecule.log_entry(note))
        print("Poly A selection step complete")
        return retained_sample

    def three_prime_bias(self, molecule, tail_length, note):
        breakage_likelihood = min([self.min_breakage_prob * (len(molecule.sequence) - tail_length), self.max_breakage_prob])
        breakage = np.random.choice([1, 0], 1, p=[breakage_likelihood, 1 - breakage_likelihood])
        if breakage:
            breakpoint_from_3p = min(np.random.geometric(self.breakpoint_prob),len(molecule.sequence) - tail_length)
            breakpoint_from_5p = len(molecule.sequence) - 1 - tail_length - breakpoint_from_3p
            molecule.truncate(breakpoint_from_5p)
            note += ' broken'
        return note

    def validate(self):
        print(f"Poly A step validating parameters")
        if self.min_retention_prob < 0 or self.min_retention_prob > 1:
            print("The minimum retention probability parameter must be between 0 and 1 inclusive", file=sys.stderr)
            return False
        if self.max_retention_prob < 0 or self.max_retention_prob > 1:
            print("The maximum retention probability parameter must be between 0 and 1 inclusive", file=sys.stderr)
            return False
        if self.max_retention_prob < self.min_retention_prob:
            print("The maximum retention probability parameter must be greater than or equal to the minimum retention probability.", file=sys.stderr)
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
            molecule = Molecule(i, sequence, 1, cigar)
            molecules.append(molecule)
            i += 1
    log_file =  "../../data/polyA_selection_step_log.txt"
    parameters = {
        "min_retention_prob": 0.02,
        "max_retention_prob": 0.98,
        "length_retention_prob": 0.04,
        "min_breakage_prob": 0.0005,
        "max_breakage_prob": 0.98,
        "breakpoint_prob": 0.001
      }
    step = PolyAStep(log_file, parameters)
    new_sample = step.execute(molecules)
    with open("../../data/polyASelectionStepTest.txt", "w+") as results:
        for molecule in new_sample:
            results.write(molecule.sequence)
