from beers_utils.molecule import Molecule
from beers.utilities.library_prep_utils import Utils


class SecondStrandPrimeStep:

    name = "Second Strand Prime Step"

    def __init__(self, log_file, parameters):
        self.history_filename = log_file
        print("Second_strand_primer_step instantiated")

    def execute(self, molecule_packet):
        """
        Starting with a basic implementation that assumes all fragments get
        primed perfectly at their 5' ends by random hexamers.
        """
        print("Second strand prime step starting")
        primer_length = 6
        primed_sample = []
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                if len(molecule.sequence) < primer_length:
                    # Too short to prime, skips out this step
                    continue
                primer_sequence = Utils.create_complement_strand(molecule.sequence[-primer_length:])
                primer_start = 1  # TODO: right start?
                primer_num = molecule.molecule_id + ".2nd_strand_primer"
                primer_molecule = Molecule(primer_num, primer_sequence, primer_start)
                molecule.bind(primer_molecule)
                primed_sample.append(molecule)
                log_file.write(molecule.log_entry("Primers bound - " + molecule.print_bound_molecules()))
        print("Second strand prime step complete")
        molecule_packet.molecules = primed_sample
        return molecule_packet

    def validate(self):
        print("Second strand prime step validating parameters")
        return True
