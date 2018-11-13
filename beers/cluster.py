from collections import namedtuple

Coordinates = namedtuple('Coordinates', ['flowcell', 'lane', 'tile', 'x', 'y'])

class Cluster:

    def __init__(self, molecule):
        self.molecule = molecule
        self.molecule_count = 1
        self.amplification_snps = []
        for index in range(len(molecule.sequence)):
            counts = { nt:0 for nt in 'GATC'}
            counts[molecule.sequence[index]] = 1
            self.amplification_snps.append(counts)

    def assign_coordinates(self, flowcell, lane, tile, x, y):
        self.coordinates = Coordinates(flowcell, lane, tile, x, y)
