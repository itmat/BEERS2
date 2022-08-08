import os
import glob
import contextlib
import collections
import pysam
from beers.cluster_packet import ClusterPacket
from beers_utils.general_utils import GeneralUtils
from beers_utils.constants import CONSTANTS
from beers.utilities.demultiplex import demultiplexer

class SAM:
    """
    The SAM object generates a SAM/BAM report for the run for a given flowcell  and for each direction in the case of
    paired end reads.  This report function is called by the controller but ONLY when all cluster packets have been
    processed.
    """

    def __init__(self, flowcell, cluster_packet_directory, sam_output_directory, sample_id, sample_barcodes):
        """
        The SAM object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).
        :param flowcell: The flowcell to which this SAM object applies.
        :param cluster_packet_directory: The location of the cluster packet files coming from the sequence pipeline.
        The assumption is the all the cluster packets are available, which is why the report generation is defered by
        the controller until the auditor determines that all cluster packets have been processed.
        :param sam_output_directory: The location where the SAM reports are filed.
        :param sample_id: ID of the sample we are making SAMs for
        :param sample_barcodes: dictionary mapping sample id to barcode as a string like f'{i5}+{i7}'
        """
        self.flowcell = flowcell
        self.cluster_packet_directory = cluster_packet_directory
        self.sam_output_directory = sam_output_directory
        self.sample_id = sample_id
        self.sample_barcodes = sample_barcodes

    def generate_report(self, reference_seqs, BAM=False, sort_by_coordinates=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the SAM files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest. The remaining clusters are sorted by their coordinates and each entry is written to
        the SAM file.

        :param reference_seqs: dictionary mapping reference names to reference sequences, used for SAM header
        :param BAM: if true, output in BAM format, else SAM (default)
        :param sort_by_coordinates: Whether to sort output by coordinates, as would typically be done with a fastq
                file from Illumina. Default: False. Setting this to True will consume considerably more memroy
                as all reads are read in at once.
        """
        sam_header = {
            "HD": { "VN": "1.0"},
            "SQ": [ {'SN': chrom_name.split()[0], 'LN': len(seq)}
                        for chrom_name, seq in reference_seqs.items() ]
        }
        chrom_list = [sq['SN'] for sq in sam_header['SQ']]

        cluster_packet_file_paths = glob.glob(f'{self.cluster_packet_directory}{os.sep}**{os.sep}*.gzip',
                                              recursive=True)
        print("Loading from", cluster_packet_file_paths)
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
            bad_barcode_file_path = {lane: os.path.join(self.sam_output_directory, f"S{self.sample_id}_unidentified_L{lane}.{'bam' if BAM else 'sam'}")
                                                    for lane in self.flowcell.lanes_to_use}
            sam_output_file_path = {lane: {barcode: os.path.join(self.sam_output_directory, f"S{sample}_L{lane}.{'bam' if BAM else 'sam'}")}
                                            for sample, barcode in self.sample_barcodes.items()
                                            for lane in  self.flowcell.lanes_to_use}
            print(f"Writing out demultiplexed sam files to:")
            for lane in sam_output_file_path.values():
                for sam in lane.values():
                    print(sam)

            bad_barcode_files = {lane: stack.enter_context(pysam.AlignmentFile(bad_barcode_file_path[lane], ('wb' if BAM else 'w'), header=sam_header))
                                                for lane in self.flowcell.lanes_to_use}
            def sam_by_barcode(lane):
                sam_files = {barcode: stack.enter_context(pysam.AlignmentFile(file_path, ('wb' if BAM else 'w'), header=sam_header))
                                for barcode, file_path in sam_output_file_path[lane].items()}
                return collections.defaultdict(
                    lambda : bad_barcode_files[lane],
                    **sam_files
                )
            sam_output_files = {lane: sam_by_barcode(lane) for lane in self.flowcell.lanes_to_use}
            demuxes = {lane: demultiplexer(output_files) for lane, output_files in sam_output_files.items()}

            #[cluster.generate_fasta_header(direction) for cluster in clusters] TODO: need this?
            for clusters in cluster_generator():
                for lane in self.flowcell.lanes_to_use:
                    lane_clusters = [cluster for cluster in clusters if cluster.lane == lane]

                    for cluster in lane_clusters:
                        paired = len(cluster.called_sequences) == 2
                        sam = demuxes[lane](cluster.called_barcode)
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
