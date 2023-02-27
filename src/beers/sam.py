import os
import glob
import contextlib
import collections
import pysam
from beers.cluster_packet import ClusterPacket
from beers_utils.general_utils import GeneralUtils
from beers_utils.constants import CONSTANTS
from beers.utilities.demultiplex import demultiplexer
from beers.flowcell import Flowcell

class SAM:
    """
    The SAM object generates a SAM/BAM report for the run for a given flowcell  and for each direction in the case of
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
        The SAM object requires the flowcell, the top level directory housing the cluster packets that have
        emerged from the sequence pipeline (they will be in the data directory under the sequence pipeline stage name),
        and the output directory for the fasta files (they will be in the data directory under the controller stage
        name).

        Parameters
        ----------
        flowcell:
            The flowcell to which this SAM object applies.
        sample_id:
            id of the sample whose SAM files we generate
        sample_barcode:
            barcodes as a string of the format '{i5}+{i7}'. Demultiplexing is done off these
        """
        self.flowcell = flowcell
        self.sample_id = sample_id
        self.sample_barcode = sample_barcode

    def generate_report(self, cluster_packet_paths, output_file_paths, bad_barcode_file_paths, reference_seqs, BAM=False, sort_by_coordinates=False):
        """
        The principal method of this object generates one or two reports depending upon whether paired end reads are
        called for.  All the information needed to create the SAM files is found in the cluster packets themselves.
        For each cluster packet, we identify whether there is one called sequence or two (paired ends).  The first
        called sequence in the list is always the forward one.  We find each cluster in the cluster packet that is
        affixed to the lanes of interest. The remaining clusters are sorted by their coordinates and each entry is written to
        the SAM file.

        Parameters
        ----------
        cluster_packet_paths:
            list of cluster packet file paths to a generate report from
        output_file_paths:
            list of sam/bam file to output to for each lane, in the same order as config.flowcell.lanes_to_use
        bad_barcode_file_paths:
            list of sam/bam files to output to for each lane, in the same order as config.flowcell.lanes_to_use
        reference_seqs:
            dictionary mapping reference names to reference sequences, used for SAM header
        BAM:
            if true, output in BAM format, else SAM (default)
        sort_by_coordinates:
            Whether to sort output by coordinates, as would typically be done with a fastq
            file from Illumina. Default: False. Setting this to True will consume considerably more memory
            as all reads are read in at once.
        """
        sam_header = {
            "HD": { "VN": "1.0"},
            "SQ": [ {'SN': chrom_name.split()[0], 'LN': len(seq)}
                        for chrom_name, seq in reference_seqs.items() ]
        }
        chrom_list = [sq['SN'] for sq in sam_header['SQ']]

        print("Loading from", cluster_packet_paths)
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
            sam_output_file_path = {lane: path for lane, path in zip(self.flowcell.lanes_to_use, output_file_paths)}
            bad_barcode_file_path = {lane: path for lane, path in zip(self.flowcell.lanes_to_use, bad_barcode_file_paths)}
            print(f"Writing out demultiplexed sam files to:")
            for sam_path in sam_output_file_path.values():
                print(sam_path)

            bad_barcode_files = {lane: stack.enter_context(pysam.AlignmentFile(bad_barcode_file_path[lane], ('wb' if BAM else 'w'), header=sam_header))
                                                for lane in self.flowcell.lanes_to_use}
            def sam_by_barcode(lane):
                sam_files = {self.sample_barcode: stack.enter_context(pysam.AlignmentFile(sam_output_file_path[lane], ('wb' if BAM else 'w'), header=sam_header))}
                return collections.defaultdict(
                    lambda : bad_barcode_files[lane],
                    **sam_files
                )
            sam_output_files = {lane: sam_by_barcode(lane) for lane in self.flowcell.lanes_to_use}
            demuxes = {lane: demultiplexer(output_files) for lane, output_files in sam_output_files.items()}

            for cluster in cluster_generator():
                paired = len(cluster.called_sequences) == 2
                sam = demuxes[cluster.lane](cluster.called_barcode)
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
