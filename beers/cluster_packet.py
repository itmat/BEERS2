from beers.sample import Sample
from beers.cluster import Cluster
import os
import resource

class ClusterPacket:

    next_cluster_packet_id = 1  # Static variable for creating increasing cluster packet id's

    def __init__(self, cluster_packet_id, sample, clusters):
        self.cluster_packet_id = cluster_packet_id
        self.sample = sample
        self.clusters = clusters

    def __str__(self):
        return f"cluster_packet_id: {self.cluster_packet_id}, sample_name: {self.sample.sample_name}, " \
               f"# of clusters: {len(self.clusters)}"

    def serialize(self, file_path):
        with open(file_path, 'wb') as obj_file:
            obj_file.write((f"#{self.cluster_packet_id}\n#{self.sample.serialize()}\n").encode())
            for cluster in self.clusters:
                obj_file.write(cluster.serialize().encode())
                # Clusters take up a variable number of lines, so we need a separator
                obj_file.write("-\n".encode())

    @staticmethod
    def deserialize(file_path):
        cluster_lines = []
        clusters = []
        with open(file_path, 'rb') as obj_file:
            for line_number, line in enumerate(obj_file):
                line = line.rstrip()
                if line_number == 0:
                    cluster_packet_id = int(line[1:].decode())
                elif line_number == 1:
                    sample = Sample.deserialize(line.decode())
                else:
                    if line.decode() == '-':
                        if cluster_lines:
                            clusters.append(Cluster.deserialize("\n".join(cluster_lines)))
                        cluster_lines = []
                    else:
                        cluster_lines.append(line.decode())
        return ClusterPacket(cluster_packet_id, sample, clusters)

    @staticmethod
    def get_serialized_cluster_packet(intermediate_directory_path, cluster_packet_filename):
        cluster_packet_file_path = os.path.join(intermediate_directory_path, cluster_packet_filename)
        cluster_packet = ClusterPacket.deserialize(cluster_packet_file_path)
        print(
            f"Intermediate loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")
        return cluster_packet