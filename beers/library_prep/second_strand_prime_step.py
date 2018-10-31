from beers.library_prep.molecule import Molecule
from beers.utilities.utils import Utils

class SecondStrandPrimeStep:

    def __init__(self, log_file, parameters):
        self.history_filename = log_file
        print("Second_strand_primer_step instantiated")

    def execute(self, sample):
        """
        Starting with a basic implementation that assumes all fragments get
        primed perfectly at their 3' ends by random hexamers.
        """
        print("Second strand prime step starting")
        primer_length = 6
        primed_sample = []
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in sample:
                primer_sequence = Utils.create_complement_strand(molecule.sequence[-primer_length:])
                primer_start = len(molecule) - primer_length
                primer_num = 1
                primer_molecule = Molecule(primer_num, primer_sequence, primer_start)
                molecule.bind(primer_molecule)
                primed_sample.append(molecule)
                log_file.write(molecule.log_entry("Primers bound - " + molecule.print_bound_molecules()))
        print("Second strand prime step complete")
        return primed_sample

    def validate(self):
        print("Second strand prime step validating parameters")
        return True
