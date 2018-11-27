import os
import pickle
from collections import namedtuple


class FastQ:

    def __init__(self, lane, cluster_packet_directory, fastq_output_directory):
        self.lane = lane
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory

    def generate_report(self):
        lane_clusters = []
        cluster_packet_filenames = [cluster_packet_filename for cluster_packet_filename in os.listdir(self.cluster_packet_directory)
                     if os.path.isfile(os.path.join(self.cluster_packet_directory, cluster_packet_filename))]
        for cluster_packet_filename in cluster_packet_filenames:
            cluster_packet_file_path = os.path.join(self.cluster_packet_directory, cluster_packet_filename)
            with open(cluster_packet_file_path, 'rb') as cluster_packet_file:
                cluster_packet = pickle.load(cluster_packet_file)
                sample_name = cluster_packet.sample.sample_name
            lane_clusters.extend([cluster for cluster in cluster_packet.clusters
                                  if cluster.coordinates.lane == self.lane])
            [cluster.generate_fasta_header(sample_name) for cluster in lane_clusters]
        sorted_clusters = sorted(lane_clusters, key=lambda lane_cluster: lane_cluster.coordinates.lane)

        fastq_output_file_path = os.path.join(self.fastq_output_directory, f"lane{self.lane}.fasta")
        with open(fastq_output_file_path, "w") as fastq_output_file:
            for cluster in sorted_clusters:
                fastq_output_file.write(cluster.header + "\n")
                fastq_output_file.write(cluster.called_sequences[0] + "\n")
                fastq_output_file.write("+\n")
                fastq_output_file.write(cluster.quality_scores[0] + "\n")


if __name__ == '__main__':
    fastq = FastQ(1, "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output/packets",
                  "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output")
    fastq.generate_report()





