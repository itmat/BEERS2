import os
import glob
import contextlib
import collections

from beers.cluster_packet import ClusterPacket
from beers.utilities.demultiplex import demultiplexer
from beers_utils.constants import CONSTANTS


class FastQ:
    """
    The FastQ object generates a FastQ report for the run for a given flowcell and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, flowcell, cluster_packet_directory, fastq_output_directory, sample_id, sample_barcodes):
        """
        The FastQ object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param flowcell: The flowcell to which this fastQ object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param fastq_output_directory: The location where the FASTQ reports are filed.  
        :param sample_id: id of the sample whose fastq files we generate
        :param sample_barcodes: dict mapping sample ids to barcodes as tuple (i5, i7). Demultiplexing is done off these
        """
        self.flowcell = flowcell
        self.cluster_packet_directory = cluster_packet_directory
        self.fastq_output_directory = fastq_output_directory
        self.sample_id = sample_id
        self.sample_barcodes = sample_barcodes

    def generate_report(self, sort_by_coordinates=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the FASTQ files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest.  A fasta header is generated internally for each of those remaining clusters
        for the given direction.  The remaining clusters are sorted by their coordinates and each entry is written to
        the FASTQ file - header, called sequence, + quality score.

        Output files are named according to barcode_S#_L#_R#.fastq specifying sample, lane and read direction numbers.

        :param sort_by_coordinates: Whether to sort output by coordinates, as would typically be done with a fastq
                file from Illumina. Default: False. Setting this to True will consume considerably more memroy
                as all reads are read in at once.
        """
        cluster_packet_file_paths = glob.glob(f'{self.cluster_packet_directory}{os.sep}**{os.sep}*.gzip',
                                              recursive=True)
        def cluster_generator():
            def inner_cluster_generator():
                for cluster_packet_file_path in cluster_packet_file_paths:
                    cluster_packet = ClusterPacket.deserialize(cluster_packet_file_path, skip_base_counts=True)
                    yield cluster_packet.clusters

            if sort_by_coordinates:
                yield from sorted(inner_cluster_generator(), key=lambda cluster: cluster.coordinates)
            else:
                yield from inner_cluster_generator()

        with contextlib.ExitStack() as stack:
            # Open all the files
            bad_barcode_file_path = {direction: {lane: os.path.join(self.fastq_output_directory, f"S{self.sample_id}_unidentified_L{lane}_R{direction}.fastq")
                                                    for lane in self.flowcell.lanes_to_use}
                                                for direction in CONSTANTS.DIRECTION_CONVENTION}
            fastq_output_file_path = {direction: {lane: {barcode: os.path.join(self.fastq_output_directory, f"S{sample}_L{lane}_R{direction}.fastq")
                                                    for sample, barcode in self.sample_barcodes.items()}
                                            for lane in  self.flowcell.lanes_to_use}
                                    for direction in CONSTANTS.DIRECTION_CONVENTION}
            print(f"Writing out demultiplexed fastq files to:")
            for direction in fastq_output_file_path.values():
                for lane in direction.values():
                    for fastq in lane.values():
                        print(fastq)

            bad_barcode_files = {direction: {lane: stack.enter_context(open(bad_barcode_file_path[direction][lane], "w"))
                                                for lane in self.flowcell.lanes_to_use}
                                            for direction in CONSTANTS.DIRECTION_CONVENTION}
            def fastq_by_barcode(direction, lane):
                fastq_files = {barcode: stack.enter_context(open(file_path,'w'))
                                    for barcode, file_path in fastq_output_file_path[direction][lane].items()}
                return collections.defaultdict(
                    lambda : bad_barcode_files[direction][lane],
                    fastq_files
                )
            fastq_output_files = {direction: {lane:  fastq_by_barcode(direction, lane)
                                                for lane in self.flowcell.lanes_to_use}
                                            for direction in CONSTANTS.DIRECTION_CONVENTION}
            demuxes = {direction: {lane: demultiplexer(fastqs) for lane, fastqs in fastq_output_files[direction].items()}
                                for direction in CONSTANTS.DIRECTION_CONVENTION}

            # Iterate through all the cluster packets
            for clusters in cluster_generator():
                for direction in CONSTANTS.DIRECTION_CONVENTION:
                    for lane in self.flowcell.lanes_to_use:
                        lane_clusters = [cluster for cluster in clusters if cluster.lane == lane]

                        for cluster in lane_clusters:
                            cluster.generate_fasta_header(direction)
                            barcode = cluster.called_barcode
                            fastq = demuxes[direction][lane](barcode)
                            fastq.write(cluster.header + "\n")
                            fastq.write(cluster.called_sequences[direction - 1] + "\n")
                            fastq.write("+\n")
                            fastq.write(cluster.quality_scores[direction - 1] + "\n")


if __name__ == '__main__':
    fastq = FastQ(1, "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output/packets",
                  "/home/crislawrence/Documents/beers_project/BEERS2.0/data/sequence/output")
    fastq.generate_report()
