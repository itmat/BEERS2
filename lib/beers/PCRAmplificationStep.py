import numpy as np
from molecule import Molecule
import sys
from utils import Utils
from timeit import default_timer as timer
from copy import copy

class PCRAmplificationStep:

    MAX_CYCLE_NUMBER = 5

    def __init__(self, log_filename, parameters):
        self.log_filename = log_filename
        self.number_cycles = parameters.get("number_cycles")
        self.sample_id_ctr = dict()

    def execute(self, sample):
        print("PCR Amplication step starting")
        self.sample_id_ctr = dict.fromkeys([molecule.molecule_id for molecule in sample], 0)
        molecules = sample
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for cycle in range(self.number_cycles):
                new_molecules = []
                for molecule in molecules:
                    new_molecule = copy(molecule)
                    self.assign_id(new_molecule, molecule.molecule_id)
                    new_molecules.append(new_molecule)
                molecules.extend(new_molecules)
            for molecule in molecules:
                log_file.write(molecule.log_entry())
        return molecules


    def assign_id(self, new_molecule, ancestor_id):
        while True:
            if ancestor_id in self.sample_id_ctr:
                new_molecule.id = str(ancestor_id) + "." + str(self.sample_id_ctr[ancestor_id])
                self.sample_id_ctr[ancestor_id] += 1
                break
            index = ancestor_id[:ancestor_id.rfind(".")]
            assert index < 0, f"cannot locate the ancestor of {new_molecule.id}"
            ancestor_id = ancestor_id[:index]

    def validate(self):
        print(f"PCR Amplification step validating parameters")
        if self.number_cycles < 0 or self.number_cycles > PCRAmplificationStep.MAX_CYCLE_NUMBER:
            print(f"The cycle number parameter must be greater than 0 and less"
                  f" than {PCRAmplificationStep.MAX_CYCLE_NUMBER}", file=sys.stderr)
            return False
        return True

if __name__ == "__main__":
    np.random.seed(100)
    molecules = Utils.convert_log_data_into_molecules("../../data/sizing_step.log")
    input_data_log_file = "../../data/pcr_amplification_step_input_data.log"
    with open(input_data_log_file, "w+") as input_data_log:
        input_data_log.write(Molecule.header)
        for rna_molecule in molecules:
            input_data_log.write(rna_molecule.log_entry())
    output_data_log_file = "../../data/pcr_amplification_step_output_data.log"
    input_parameters = {
        "number_cycles": 2
    }
    step = PCRAmplificationStep(output_data_log_file, input_parameters)
    start = timer()
    step.execute(molecules)
    end = timer()
    print(f"PCRAmplification Step: {end - start} for {len(molecules)} molecules input.")
