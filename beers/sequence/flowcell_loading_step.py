import pysam
import re
from beers.cluster import Coordinates


class FlowcellLoadingStep:

    coords_match_pattern = "^.*:(\w+):(\d+)\:(\d+)\:(\d+)\:(\d+)$"

    def __init__(self, cluster_packets, parameters):
        self.cluster_packets = cluster_packets
        self.input_file_paths = [cluster_packet.sample.input_file_path for cluster_packet in cluster_packets]
        self.coordinate_generator = self.generate_coordinates()

    def generate_coordinates(self):
        for input_file_path in self.input_file_paths:
            input_file = pysam.AlignmentFile(input_file_path, "rb")
            for line in input_file.fetch():
                if line.is_unmapped or not line.is_read1 or line.get_tag(tag="NH") != 1:
                    continue
                coords_match = re.match(FlowcellLoadingStep.coords_match_pattern, line.qname)
                if coords_match:
                    (flowcell, lane, tile, x, y) = coords_match.groups()
                    coordinates = Coordinates(flowcell, lane, tile, x, y)
                    yield(coordinates)

    def execute(self):
        for cluster_packet in self.cluster_packets:
            for cluster in cluster_packet.clusters:
                cluster.assign_coordinates(next(self.coordinate_generator))
        return self.cluster_packets
