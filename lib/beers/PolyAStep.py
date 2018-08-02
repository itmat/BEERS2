from molecule import Molecule
import numpy as np


class PolyAStep:

    def __init__(self):
        self.history_filename = "polyA_selection_history.txt"
        print("Poly A selection step instantiated")

    def execute(self, sample):
        print("Poly A selection step starting")
        retained_sample = []
        with open(self.history_filename, "w+") as history_file:
            for molecule in sample:
                tail_length = molecule.poly_a_tail_length()
                retention_odds = min(0.02 + 0.01 * tail_length, 0.99)
                retained = np.random.choice([1, 0], 1, p=[retention_odds, 1 - retention_odds])[0]
                note = ''
                if retained:
                    retained_sample.append(molecule)
                    note += 'retained'
                    breakage_frequency = min(0.01 * len(molecule.sequence) - tail_length, 0.99)
                    breakage = np.random.choice([1,0], 1, p=[breakage_frequency, 1 - breakage_frequency])
                    if breakage:
                        break_distribution = np.random.geometric(0.01, len(molecule.sequence) - tail_length)
                        break_distribution
                else:
                    note += 'removed'
                #print(tail_length, retention_odds, note)
                history_file.write(molecule.log_entry() + "," + note + "\n")
        print("Poly A selection step complete")
        return retained_sample

    def validate(self, **kwargs):
        print(f"Poly A step validating input kw args")
        return True


if __name__ == "__main__":
    molecules = []
    i = 0
    with open("molecules.txt", 'r') as sample_file:
        sequences = sample_file.readlines()
        for sequence in sequences:
            cigar = str(len(sequence)) + 'M'
            molecule = Molecule(i, sequence, '1', cigar)
            molecules.append(molecule)
            i += 1
    step = PolyAStep()
    new_sample = step.execute(molecules)
    with open("polyASelectionStepTest.txt", "w+") as results:
        for molecule in new_sample:
            results.write(molecule.sequence)
