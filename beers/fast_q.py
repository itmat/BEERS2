import os
import glob
from beers.cluster_packet import ClusterPacket


class FastQ:

    def __init__(self, lane, cluster_packet_directory, fastq_output_directory, forward=True):
        self.lane = lane
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory

    def generate_report(self):
        lane_clusters = []
        cluster_packet_file_paths = glob.glob(f'{self.cluster_packet_directory}{os.sep}**{os.sep}*.gzip',
                                              recursive=True)
        for cluster_packet_file_path in cluster_packet_file_paths:
            cluster_packet = ClusterPacket.deserialize(cluster_packet_file_path)
            sample_name = cluster_packet.sample.sample_name
            lane_clusters += [cluster for cluster in cluster_packet.clusters if cluster.lane == self.lane]
            [cluster.generate_fasta_header("f") for cluster in lane_clusters]
        sorted_clusters = sorted(lane_clusters, key=lambda lane_cluster: lane_cluster.lane)

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





