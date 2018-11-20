class MoleculePacket:

    next_molecule_packet_id = 1  # Static variable for creating increasing molecule packet id's

    def __init__(self, molecule_packet_id, sample, molecules):
        self.molecule_packet_id = molecule_packet_id
        self.sample = sample
        self.molecules = molecules
