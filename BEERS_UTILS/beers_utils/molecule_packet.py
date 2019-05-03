from beers_utils.sample import Sample
from beers_utils.molecule import Molecule
import os
import resource

class MoleculePacket:

    next_molecule_packet_id = 0  # Static variable for creating increasing molecule packet id's

    def __init__(self, molecule_packet_id, sample, molecules):
        self.molecule_packet_id = molecule_packet_id
        self.sample = sample
        self.molecules = molecules

    def serialize(self, file_path):
        with open(file_path, 'wb') as obj_file:
            obj_file.write((f"#{self.molecule_packet_id}\n#{self.sample.serialize()}\n").encode(encoding="ascii"))
            for molecule in self.molecules:
                obj_file.write((molecule.serialize() + "\n").encode(encoding="ascii"))

    @staticmethod
    def deserialize(file_path):
        molecules = []
        with open(file_path, 'rb') as obj_file:
            for line_number, line in enumerate(obj_file):
                line = line.rstrip()
                if line_number == 0:
                    molecule_packet_id = int(line[1:].decode(encoding="ascii"))
                elif line_number == 1:
                    sample = Sample.deserialize(line.decode(encoding="ascii"))
                else:
                    molecules.append(Molecule.deserialize(line.decode(encoding="ascii")))
        return MoleculePacket(molecule_packet_id, sample, molecules)

    @staticmethod
    def get_serialized_molecule_packet(input_directory_path, molecule_packet_filename):
        molecule_packet_file_path = os.path.join(input_directory_path, molecule_packet_filename)
        molecule_packet = MoleculePacket.deserialize(molecule_packet_file_path)
        print(
            f"Input loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")
        return molecule_packet
