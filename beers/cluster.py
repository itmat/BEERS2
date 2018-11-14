from collections import namedtuple
import numpy as np

BaseCounts = namedtuple('BaseCounts', ["G", "A", "T", "C"])
Coordinates = namedtuple('Coordinates', ['flowcell', 'lane', 'tile', 'x', 'y'])

class Cluster:

    next_cluster_id = 1  # Static variable for creating increasing cluster id's

    def __init__(self, cluster_id, molecule):
        self.cluster_id = cluster_id
        self.molecule = molecule
        self.molecule_count = 1
        counts = {nt: np.array([1 if base == nt else 0 for base in molecule.sequence]) for nt in "ACGT"}
        self.base_counts = BaseCounts(counts['G'],counts['A'], counts['T'], counts['C'])

    def assign_coordinates(self, flowcell, lane, tile, x, y):
        self.coordinates = Coordinates(flowcell, lane, tile, x, y)
