import os
import io
import glob
import contextlib
import collections

from beers.cluster_packet import ClusterPacket
from beers.utilities.demultiplex import demultiplexer
from beers_utils.constants import CONSTANTS
from beers.flowcell import Flowcell


class FastQ:
    """
    The FastQ object generates a FastQ report for the run for a given flowcell and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(
            self,
            flowcell: Flowcell,
            sample_id: str,
            sample_barcode: str,
        ):
        """
        The FastQ object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).

        Parameters
        ----------
        flowcell:
            The flowcell to which this fastQ object applies.
        sample_id:
            id of the sample whose fastq files we generate
        sample_barcode:
            barcodes as a string of the format '{i5}+{i7}'. Demultiplexing is done off these
        """
        self.flowcell = flowcell
        self.sample_id = sample_id
        self.sample_barcode = sample_barcode

    def generate_report(self, cluster_packet_paths, output_file_paths, bad_barcode_file_paths, sort_by_coordinates=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the FASTQ files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest.  A fasta header is generated internally for each of those remaining clusters
        for the given direction.  The remaining clusters are sorted by their coordinates and each entry is written to
        the FASTQ file - header, called sequence, + quality score.

        Output files are named according to barcode_S#_L#_R#.fastq specifying sample, lane and read direction numbers.


        Parameters
        ----------
        cluster_packet_paths:
            list of cluster packet file paths to a generate report from
        output_file_paths:
            list of two (one for R1/R2 read directions) lists of sam/bam file to output to for each lane, in the same order as config.flowcell.lanes_to_use
        bad_barcode_file_paths:
            list of two (one for R1/R2 read directions) list of sam/bam files to output to for each lane, in the same order as config.flowcell.lanes_to_use
        sort_by_coordinates:
                Whether to sort output by coordinates, as would typically be done with a fastq
                file from Illumina. Default: False. Setting this to True will consume considerably more memory
                as all reads are read in at once.
        """
        def cluster_generator():
            def inner_cluster_generator():
                for cluster_packet_path in cluster_packet_paths:
                    cluster_packet = ClusterPacket.deserialize(cluster_packet_path, skip_base_counts=True)
                    yield from cluster_packet.clusters

            if sort_by_coordinates:
                yield from sorted(inner_cluster_generator(), key=lambda cluster: cluster.coordinates)
            else:
                yield from inner_cluster_generator()

        with contextlib.ExitStack() as stack:
            # Open all the files
            fastq_output_file_path = {direction: {lane: path for lane, path in zip(self.flowcell.lanes_to_use, out_paths)}
                                        for direction, out_paths in zip(CONSTANTS.DIRECTION_CONVENTION, output_file_paths)}
            bad_barcode_file_path = {direction: {lane: path for lane, path in zip(self.flowcell.lanes_to_use, bad_paths)}
                                        for direction, bad_paths in zip(CONSTANTS.DIRECTION_CONVENTION, bad_barcode_file_paths)}
            print(f"Writing out demultiplexed fastq files to:")
            for direction in fastq_output_file_path.values():
                for fastq in direction.values():
                    print(fastq)

            bad_barcode_files = {direction: {lane: stack.enter_context(open(bad_barcode_file_path[direction][lane], "w"))
                                                for lane in self.flowcell.lanes_to_use}
                                            for direction in CONSTANTS.DIRECTION_CONVENTION}
            def fastq_by_barcode(direction, lane):
                fastq_files = {self.sample_barcode: stack.enter_context(open(fastq_output_file_path[direction][lane],'w'))}
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
            for cluster in cluster_generator():
                for direction_num in CONSTANTS.DIRECTION_CONVENTION:
                    cluster.generate_fasta_header(direction_num)
                    barcode = cluster.called_barcode
                    fastq_file = demuxes[direction_num][cluster.lane](barcode)
                    fastq_file.write(cluster.header + "\n")
                    fastq_file.write(cluster.called_sequences[direction_num - 1] + "\n")
                    fastq_file.write("+\n")
                    fastq_file.write(cluster.quality_scores[direction_num - 1] + "\n")
