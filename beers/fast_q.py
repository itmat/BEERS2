import os
import glob
import contextlib
import collections

from beers.cluster_packet import ClusterPacket
from beers_utils.constants import CONSTANTS


class FastQ:
    """
    The FastQ object generates a FastQ report for the run for a given flowcell and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, flowcell, cluster_packet_directory, fastq_output_directory, sample_barcodes):
        """
        The FastQ object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param flowcell: The flowcell to which this fastQ object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param fastq_output_directory: The location where the FASTQ reports are filed.  Note that no organization into
        subdirectories is needed here since compartively few reports are generated (at most 2 per flowcell lane).
        :param sample_barcodes: dict mapping sample ids to barcodes as tuple (i5, i7). Demultiplexing is done off these
        """
        self.flowcell = flowcell
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory
        self.sample_barcodes = sample_barcodes

    def generate_report(self):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the FASTQ files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest.  A fasta header is generated internally for each of those remaining clusters
        for the given direction.  The remaining clusters are sorted by their coordinates and each entry is written to
        the FASTQ file - header, called sequence, + quality score.

        Output files are named according to barcode_S#_L#_R#.fastq specifying sample, lane and read direction numbers.
        """
        cluster_packet_file_paths = glob.glob(f'{self.cluster_packet_directory}{os.sep}**{os.sep}*.gzip',
                                              recursive=True)
        clusters = []
        max_direction_num = min(CONSTANTS.DIRECTION_CONVENTION)
        for cluster_packet_file_path in cluster_packet_file_paths:
            cluster_packet = ClusterPacket.deserialize(cluster_packet_file_path)
            max_direction_num = max(max_direction_num, len(cluster_packet.clusters[0].called_sequences))
            clusters += cluster_packet.clusters

        for direction in CONSTANTS.DIRECTION_CONVENTION:
            if direction > max_direction_num:
                break # Processed all available read directions, nothing else to do

            [cluster.generate_fasta_header(direction) for cluster in clusters]
            for lane in self.flowcell.lanes_to_use:
                lane_clusters = [cluster for cluster in clusters if cluster.lane == lane]
                sorted_clusters = sorted(lane_clusters, key=lambda cluster: cluster.coordinates)

                fastq_output_file_path = {barcode: os.path.join(self.fastq_output_directory, f"S{sample}_L{lane}_R{direction}.fastq")
                                                for sample, barcode in self.sample_barcodes.items()}
                bad_barcode_file_path = os.path.join(self.fastq_output_directory, f"unidentified_L{lane}_R{direction}.fastq")
                with contextlib.ExitStack() as stack:
                    print(f"Writing out demultiplexed fastq files to: {list(fastq_output_file_path.values())}")
                    bad_barcode_file = stack.enter_context(open(bad_barcode_file_path, "w")) # Write to this file if no matching barcode
                    files = collections.defaultdict(lambda : bad_barcode_file)
                    files.update({barcode: stack.enter_context(open(fastq, "w")) for barcode, fastq in fastq_output_file_path.items()})
                    for cluster in sorted_clusters:
                        barcode = cluster.called_barcode
                        fastq = files[barcode]
                        fastq.write(cluster.header + "\n")
                        fastq.write(cluster.called_sequences[direction - 1] + "\n")
                        fastq.write("+\n")
                        fastq.write(cluster.quality_scores[direction - 1] + "\n")


if __name__ == '__main__':
    fastq = FastQ(1, "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output/packets",
                  "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output")
    fastq.generate_report()
