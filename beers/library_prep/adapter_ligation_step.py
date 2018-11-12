from beers.molecule import Molecule
import numpy as np
import pickle
from timeit import default_timer as timer


class AdapterLigationStep:

    def __init__(self, step_log_file_path, parameters):
        self.log_filename = step_log_file_path
        self.parameters = parameters
        print("Adapter Ligation step instantiated")

    def execute(self, molecule_packet):
        print("Adapter ligation step starting")
        sample = molecule_packet.sample
        with open(self.log_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                sequence = molecule.sequence
                cigar = f"{molecule.cigar or len(sequence)}M"
                molecule.sequence = sample.adapter_5_prime + sequence + sample.adapter_3_prime
                molecule.cigar = f"{len(sample.adapter_5_prime)}S{cigar}{len(sample.adapter_3_prime)}S"
                log_file.write(molecule.log_entry())
        return molecule_packet

    def validate(self):
        print(f"Adapter ligation step validating parameters")
        return True

if __name__ == "__main__":
    np.random.seed(100)
    with open("../../data/tests/molecule_packet.pickle", 'rb') as molecule_packet_file:
        molecule_packet = pickle.load(molecule_packet_file)
    input_data_log_file = "../../data/tests/adapter_ligation_step_input_data.log"
    with open(input_data_log_file, "w+") as input_data_log:
        input_data_log.write(Molecule.header)
        for rna_molecule in molecule_packet.molecules:
            input_data_log.write(rna_molecule.log_entry())
    step_log_file_path = "../../data/tests/adapter_ligation_step_output_data.log"
    input_parameters = {}
    step = AdapterLigationStep(step_log_file_path, input_parameters)
    start = timer()
    step.execute(molecule_packet)
    end = timer()
    print(f"Adapter Ligation Step: {end - start} for {len(molecule_packet.molecules)} molecules.")