from beers_utils.molecule import Molecule
from beers.utilities.library_prep_utils import Utils


class FirstStrandSynthesisStep:

    name = "First Strand Synthsis Step"

    def __init__(self, log_file, parameters):
        self.history_filename = log_file
        print("First Strand cDNA Synthesis Step instantiated.")

    def execute(self, molecule_packet):
        """



        """
        print("First strand cDNA synthesis step starting...")
        cdna_sample = []
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                cdna_seq = Utils.create_complement_strand(molecule.sequence)

                for primer in molecule.bound_molecules:
                    start = len(molecule.sequence) - (primer.start + len(primer))
                    end = start + len(primer)
                    cdna_seq = cdna_seq[:start] + primer.sequence + cdna_seq[end:]
                    #cdna_seq[start:end]  = primer.sequence

                first_primer = molecule.bound_molecules[0]
                cdna_start = len(molecule.sequence) - (first_primer.start + len(first_primer) )
                cdna_molecule = Molecule(molecule.molecule_id + '.cdna' ,cdna_seq[cdna_start:],0)
                cdna_sample.append(cdna_molecule)
                log_file.write(cdna_molecule.log_entry())

            print("First strand cDNA synthesis step complete.")
            molecule_packet.molecules = cdna_sample
            return molecule_packet



    def validate(self):
        return True
