import pysam
import numpy as np
import math
import re
from beers.molecule_packet import MoleculePacket
from beers.cluster import Cluster
from beers.cluster_packet import ClusterPacket
from beers.cluster import Coordinates

class FlowcellLoader:

    coords_match_pattern = "^.*:(\w+):(\d+)\:(\d+)\:(\d+)\:(\d+)$"

    def __init__(self, molecule_packet, parameters):
        self.molecule_packet = molecule_packet
        self.input_file_path = molecule_packet.sample.input_file_path
        self.parameters = parameters
        self.flowcell_retention = self.parameters["flowcell_retention_percentage"]/100
        self.coordinate_generator = self.generate_coordinates()

    def identify_retained_molecules(self):
        number_samples_to_draw = math.floor(self.flowcell_retention * len(self.molecule_packet.molecules))
        return np.random.choice(self.molecule_packet.molecules, size=number_samples_to_draw, replace=False)

    @staticmethod
    def convert_molecule_pkt_to_cluster_pkt(molecule_packet):
        clusters = []
        for molecule in molecule_packet.molecules:
            cluster_id = Cluster.next_cluster_id
            clusters.append(Cluster(cluster_id, molecule))
            Cluster.next_cluster_id += 1
        return ClusterPacket(molecule_packet.sample, clusters)

    def generate_coordinates(self):
        input_file = pysam.AlignmentFile(self.input_file_path, "rb")
        for line in input_file.fetch():
            if line.is_unmapped or not line.is_read1 or line.get_tag(tag="NH") != 1:
                continue
            coords_match = re.match(FlowcellLoader.coords_match_pattern, line.qname)
            if coords_match:
                (flowcell, lane, tile, x, y) = coords_match.groups()
                coordinates = Coordinates(flowcell, lane, tile, x, y)
                yield(coordinates)

    def load_flowcell(self):
        retained_molecules = self.identify_retained_molecules()
        retained_molecule_packet = MoleculePacket(MoleculePacket.next_molecule_packet_id,
                                                  self.molecule_packet.sample, retained_molecules)
        MoleculePacket.next_molecule_packet_id += 1
        cluster_packet = self.convert_molecule_pkt_to_cluster_pkt(retained_molecule_packet)
        for cluster in cluster_packet.clusters:
            cluster.assign_coordinates(next(self.coordinate_generator))
        return cluster_packet