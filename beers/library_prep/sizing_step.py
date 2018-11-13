import sys
from timeit import default_timer as timer
from scipy.stats import norm
import numpy as np
import pylab as pl
import pickle
from beers.utilities.library_prep_utils import Utils

from beers.molecule import Molecule

class SizingStep:

    def __init__(self, step_log_file_path, parameters):
        self.idealized = False
        self.log_filename = step_log_file_path
        self.mean_length = parameters.get("mean_length")
        self.sd_length = parameters.get("sd_length")
        if not self.mean_length or not self.sd_length:
            self.idealized = True
            self.min_length = parameters.get("min_length")
            self.max_length = parameters.get("max_length")
        else:
            self.dist = norm(self.mean_length, self.sd_length)
            mean_prob = self.dist.pdf(self.mean_length)
            self.normalize_coefficient = 1 / mean_prob
            self.lower_breakpoint = self.mean_length + 0.25 * self.sd_length
            self.upper_breakpoint = self.mean_length + 6 * self.sd_length
        print("Sizing step instantiated")

    def execute(self, molecule_packet):
        print("Sizing step starting")
        retained_molecules = []
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                seq_length = len(molecule.sequence)
                if self.idealized:
                    retained = self.min_length <= seq_length <= self.max_length
                else:
                    retention_odds = self.dist_function(seq_length)
                    retained = (np.random.random() < retention_odds)
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

    def dist_function(self, x):
        y = self.normalize_coefficient * self.dist.pdf(x)
        if self.lower_breakpoint < x <= self.upper_breakpoint:
                y += max(0, ((x - self.lower_breakpoint) / self.lower_breakpoint) * (1.0 - x/self.upper_breakpoint))
        return y

    def display_dist_function(self):
        x_values = np.linspace(0, self.mean_length + 6*self.sd_length, self.mean_length + 6*self.sd_length)
        y_values = [self.dist_function(x) for x in x_values]
        pl.plot(x_values, y_values)
        pl.show()


    def validate(self):
        print(f"Sizing step validating parameters")
        return True

if __name__ == "__main__":
    np.random.seed(100)

    # Getting original molecule packet (to preserve original sample metadata in case it is needed)
    with open("../../data/tests/molecule_packet.pickle", 'rb') as molecule_packet_file:
        molecule_packet = pickle.load(molecule_packet_file)

    # Taking advantage of an existing log file to grab molecules.
    molecule_packet.molecules = \
        Utils.convert_log_data_into_molecules("../../data/tests/polya_step_output_data.log")

    # Copying these molecules into a separate log file
    input_data_log_file = "../../data/tests/sizing_step_input_data.log"
    with open(input_data_log_file, "w+") as input_data_log:
        input_data_log.write(Molecule.header)
        for rna_molecule in molecule_packet.molecules:
            input_data_log.write(rna_molecule.log_entry())

    # Selecting step log file and parameter info and using both to instantiate a step
    # object (not bothering with validation)
    output_data_log_file = "../../data/tests/sizing_step_output_data.log"
    input_parameters = {
        "min_length": 100,
        "max_length": 400
    }
    step = SizingStep(output_data_log_file, input_parameters)
    #step.display_dist_function()

    # Executing the step and noting the time taken.
    start = timer()
    step.execute(molecule_packet)
    end = timer()
    print(f"Sizing Step: {end - start} for {len(molecule_packet.molecules)} molecules.")
