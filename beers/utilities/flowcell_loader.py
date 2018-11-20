import pysam
import numpy as np
import math
import re
import os
import gzip
import itertools
import sys
import json
from timeit import default_timer as timer

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

    @staticmethod
    def read_fastq(configuration):
        min_coords = {"lane": 10000, "tile": 10000, "x": 10000, "y": 10000}
        max_coords = {"lane": 0, "tile": 0, "x": 0, "y": 0}
        input_directory_path = configuration['input']['directory_path']
        for fastq_filename in configuration['input']['fastq_filenames']:
            fastq_file_path = os.path.join(input_directory_path, fastq_filename)
            with gzip.open(fastq_file_path, 'rb') as fastq_file:
                for counter, (byte_header, content, toss, scores) in enumerate(itertools.zip_longest(*[fastq_file] * 4)):
                    if counter%1000000 == 0:
                        print(f"{counter} reads")
                    header = byte_header.decode()
                    sequence_identifier = header.split(" ")[0]
                    coords_match = re.match(FlowcellLoader.coords_match_pattern, sequence_identifier)
                    if coords_match:
                        (flowcell, lane, tile, x, y) = coords_match.groups()
                        lane, tile, x, y = int(lane), int(tile), int(x), int(y)
                        min_coords['x'] = min(x, min_coords['x'])
                        max_coords['x'] = max(x, max_coords['x'])
                        min_coords['y'] = min(y, min_coords['y'])
                        max_coords['y'] = max(y, max_coords['y'])
                        min_coords['tile'] = min(tile, min_coords['tile'])
                        max_coords['tile'] = max(tile, max_coords['tile'])
                        min_coords['lane'] = min(lane, min_coords['lane'])
                        max_coords['lane'] = max(lane, max_coords['lane'])
        return min_coords, max_coords

    @staticmethod
    def create_coordinate_list():
        x = list(range(1012, 32816))
        np.random.shuffle(x)
        print(x[:25])
        y = list(range(998, 49247))
        np.random.shuffle(y)
        print(y[:25])
        tile = list(range(1101, 2228))
        np.random.shuffle(tile)
        print(tile[:25])
        lanes = [1, 3, 7]
        np.random.shuffle(lanes)
        print(lanes)
        values = itertools.product(lanes, tile, x, y)
        for value in values:
            print(value)
            input()
        #coord_list = [value for value in values]
        #print(sys.getsizeof(coord_list))



if __name__ == '__main__':
    start = timer()
    #with open("../../config/config.json", "r+") as configuration_file:
    #    configuration = json.load(configuration_file)
    #min_coords, max_coords = FlowcellLoader.read_fastq(configuration["controller"])
    #print(min_coords)
    #print(max_coords)
    FlowcellLoader.create_coordinate_list()
    end = timer()
    print(f"FastQ query: {end - start}")