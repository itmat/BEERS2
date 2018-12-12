import os
import glob
from beers.cluster_packet import ClusterPacket
from beers.constants import CONSTANTS


class FastQ:

    def __init__(self, lane, cluster_packet_directory, fastq_output_directory, forward=True):
        self.lane = lane
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory

    def generate_report(self):
        cluster_packet_file_paths = glob.glob(f'{self.cluster_packet_directory}{os.sep}**{os.sep}*.gzip',
                                                 recursive=True)
        for direction in CONSTANTS.DIRECTION_CONVENTION:
            lane_clusters = []
            abort = False
            for cluster_packet_file_path in cluster_packet_file_paths:
                cluster_packet = ClusterPacket.deserialize(cluster_packet_file_path)
                if len(cluster_packet.clusters[0].called_sequences) < direction:
                    abort = True
                    break
                sample_name = cluster_packet.sample.sample_name
                lane_clusters += [cluster for cluster in cluster_packet.clusters if cluster.lane == self.lane]
                [cluster.generate_fasta_header(direction) for cluster in lane_clusters]
            sorted_clusters = sorted(lane_clusters, key=lambda lane_cluster: lane_cluster.lane)

            fastq_output_file_path = os.path.join(self.fastq_output_directory,
                                                  f"lane{self.lane}_{direction}.fasta")
            if abort:
                break
            with open(fastq_output_file_path, "w") as fastq_output_file:
                for cluster in sorted_clusters:
                    fastq_output_file.write(cluster.header + "\n")
                    fastq_output_file.write(cluster.called_sequences[direction - 1] + "\n")
                    fastq_output_file.write("+\n")
                    fastq_output_file.write(cluster.quality_scores[direction - 1] + "\n")


if __name__ == '__main__':
    fastq = FastQ(1, "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output/packets",
                  "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output")
    fastq.generate_report()





