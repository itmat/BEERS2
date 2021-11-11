import os
import glob
import pysam
from beers.cluster_packet import ClusterPacket
from beers_utils.general_utils import GeneralUtils
from beers_utils.constants import CONSTANTS


class SAM:
    """
    The SAM object generates a SAM/BAM report for the run for a given lane and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, lane, cluster_packet_directory, sam_output_directory):
        """
        The SAM object requires the flowcell lane, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param lane: The flowcell to which this SAM object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param sam_output_directory: The location where the SAM reports are filed.  Note that no organization into
        subdirectories is needed here since compartively few reports are generated.
        """
        self.lane = lane
        self.cluster_packet_directory = cluster_packet_directory
        self.sam_output_directory = sam_output_directory

    def generate_report(self, reference_seqs, BAM=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the SAM files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lane of interest. The remaining clusters are sorted by their coordinates and each entry is written to
        the SAM file.

        :param reference_seqs: dictionary mapping reference names to reference sequences, used for SAM header
        :param BAM: if true, output in BAM format, else SAM (default)
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

            sam_output_file_path = os.path.join(self.sam_output_directory,
                                                  f"lane{self.lane}.{'bam' if BAM else 'sam'}")
            if abort:
                break

            print(f"Writing out SAM file to {sam_output_file_path}")
            sam_header = {
                "HD": { "VN": "1.0"},
                "SQ": [ {'SN': chrom_name.split()[0], 'LN': len(seq)}
                            for chrom_name, seq in reference_seqs.items() ]
            }

            chrom_list = [sq['SN'] for sq in sam_header['SQ']]

            with pysam.AlignmentFile(sam_output_file_path, ('wb' if BAM else 'w'), header=sam_header) as sam:
                for cluster in sorted_clusters:
                    paired = len(cluster.called_sequences) == 2
                    for direction, (seq, qual, start, cigar) in enumerate(zip(cluster.called_sequences, cluster.quality_scores, cluster.read_starts, cluster.read_cigars)):
                        a = pysam.AlignedSegment()
                        a.query_name = cluster.encode_sequence_identifier()
                        # TODO: are these the right flags?
                        rev_strand = ((cluster.molecule.source_strand == '-' and direction == 0) or (cluster.molecule.source_strand == '+' and direction == 1))
                        a.flag = (0x01*paired) + 0x02 + (0x40 if (direction == 0) else 0x80) + (0x10 if rev_strand else 0x20)
                        a.query_sequence = seq if not rev_strand else GeneralUtils.create_complement_strand(seq)
                        a.reference_id = chrom_list.index(cluster.molecule.source_chrom)
                        a.reference_start = start - 1 # pysam uses 0-based index, we use 1-based
                        a.mapping_quality = 255
                        a.cigarstring = cigar
                        a.query_qualities = pysam.qualitystring_to_array(qual)
                        sam.write(a)
