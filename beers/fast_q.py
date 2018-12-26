import os
import glob
from beers.cluster_packet import ClusterPacket
from beers.constants import CONSTANTS


class FastQ:
    """
    The FastQ object generates a FastQ report for the run for a given lane and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, lane, cluster_packet_directory, fastq_output_directory):
        """
        The FastQ object requires the flowcell lane, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param lane: The flowcell to which this fastQ object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param fastq_output_directory: The location where the FASTQ reports are filed.  Note that no organization into
        subdirectories is needed here since compartively few reports are generated (at most 2 per flowcell lane).
        """
        self.lane = lane
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory

    def generate_report(self):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the FASTQ files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lane of interest.  A fasta header is generated internally for each of those remaining clusters
        for the given direction.  The remaining clusters are sorted by their coordinates and each entry is written to
        the FASTQ file - header, called sequence, + quality score.
        """
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
                lane_clusters += [cluster for cluster in cluster_packet.clusters if cluster.lane == self.lane]
                [cluster.generate_fasta_header(direction) for cluster in lane_clusters]
            sorted_clusters = sorted(lane_clusters, key=lambda lane_cluster: lane_cluster.coordinates)

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
