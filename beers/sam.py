import os
import glob
import contextlib
import collections
import pysam
from beers.cluster_packet import ClusterPacket
from beers_utils.general_utils import GeneralUtils
from beers_utils.constants import CONSTANTS


class SAM:
    """
    The SAM object generates a SAM/BAM report for the run for a given flowcell  and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, flowcell, cluster_packet_directory, sam_output_directory, sample_barcodes):
        """
        The SAM object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param flowcell: The flowcell to which this SAM object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param sam_output_directory: The location where the SAM reports are filed.  Note that no organization into
        subdirectories is needed here since compartively few reports are generated.
        :param sample_barcodes: dict mapping sample ids to barcodes as tuple (i5, i7). Demultiplexing is done off these
        """
        self.flowcell = flowcell
        self.cluster_packet_directory = cluster_packet_directory
        self.sam_output_directory = sam_output_directory
        self.sample_barcodes = sample_barcodes

    def generate_report(self, reference_seqs, BAM=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the SAM files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest. The remaining clusters are sorted by their coordinates and each entry is written to
        the SAM file.

        :param reference_seqs: dictionary mapping reference names to reference sequences, used for SAM header
        :param BAM: if true, output in BAM format, else SAM (default)
        """
        sam_header = {
            "HD": { "VN": "1.0"},
            "SQ": [ {'SN': chrom_name.split()[0], 'LN': len(seq)}
                        for chrom_name, seq in reference_seqs.items() ]
        }
        chrom_list = [sq['SN'] for sq in sam_header['SQ']]

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

                sam_output_file_paths = {barcode: os.path.join(self.sam_output_directory,
                                                      f"S{sample}_L{lane}.{'bam' if BAM else 'sam'}")
                                                for sample, barcode in self.sample_barcodes.items()}
                bad_barcode_file_path = os.path.join(self.sam_output_directory, f"unidentified_L{lane}.{'bam' if BAM else 'sam'}")


                with contextlib.ExitStack() as stack:
                    print(f"Writing out demultiplexed alignment files to: {list(sam_output_file_paths.values())}")
                    # sam/bam file with reads that were not demultiplexed
                    bad_barcode_file = stack.enter_context(pysam.AlignmentFile(bad_barcode_file_path, ('wb' if BAM else 'w'), header=sam_header))
                    files = collections.defaultdict(lambda : bad_barcode_file)
                    files.update({barcode: stack.enter_context(pysam.AlignmentFile(file_path, ('wb' if BAM else 'w'), header=sam_header))
                                for barcode, file_path in sam_output_file_paths.items()})
                    for cluster in sorted_clusters:
                        paired = len(cluster.called_sequences) == 2
                        sam = files[cluster.called_barcode]
                        for direction, (seq, qual, start, cigar) in enumerate(zip(cluster.called_sequences, cluster.quality_scores, cluster.read_starts, cluster.read_cigars)):
                            a = pysam.AlignedSegment()
                            a.query_name = cluster.encode_sequence_identifier()
                            rev_strand = ((cluster.molecule.source_strand == '-' and direction == 0) or (cluster.molecule.source_strand == '+' and direction == 1))
                            a.flag = (0x01*paired) + 0x02 + (0x40 if (direction == 0) else 0x80) + (0x10 if rev_strand else 0x20)
                            a.query_sequence = seq if not rev_strand else GeneralUtils.create_complement_strand(seq)
                            a.reference_id = chrom_list.index(cluster.molecule.source_chrom)
                            a.reference_start = start - 1 # pysam uses 0-based index, we use 1-based
                            a.mapping_quality = 255
                            a.cigarstring = cigar
                            a.query_qualities = pysam.qualitystring_to_array(qual)
                            sam.write(a)
