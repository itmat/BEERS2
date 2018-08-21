import sys
from molecule import Molecule
from scipy.stats import binom
import numpy as np

class SizingStep:

    def __init__(self, log_filename, parameters):
        self.log_filename = log_filename
        self.size_mean = parameters.get("size_mean")
        self.size_sd = parameters.get("size_sd")
        var = (self.size_sd)^2
        self.p = 1 - var/self.size_mean
        self.n = self.size_mean/self.p
        print("Sizing step instantiated")

    def execute(self, sample):
        print("Sizing step starting")
        retained_sample = []
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in sample:
                seq_length = len(molecule.sequence)
                retention_odds = binom.pmf(seq_length, self.n, self.p)
                retained = np.random.choice([1, 0], 1, p=[retention_odds, 1 - retention_odds])[0]
                note = ''
                if retained:
                    retained_sample.append(molecule)
                    note += 'retained'
                else:
                    note += 'removed'
                log_file.write(molecule.log_entry(note))
        print("Sizing step complete")
        return retained_sample


    def validate(self):
        print(f"Sizing step validating parameters")
        if self.size_mean < 0 or self.size_mean > 10000:
            print("The mean size for the molecules processed in this step must be between 0 and 10000 exclusive", file=sys.stderr)
            return False
        if self.size_sd < 0 or self.size_sd > self.size_mean:
            print("The sizing standard deviation for the moleules processed in thsi step must be betwwen 0 and the given mean, exclusive", file=sys.stderr)
            return False
        return True