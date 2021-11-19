from collections import namedtuple
import numpy as np
from io import StringIO
import math
from contextlib import closing
from beers_utils.general_utils import GeneralUtils
from beers_utils.molecule import Molecule
from beers_utils.constants import CONSTANTS
import beers_utils.cigar


# BaseCounts - tuple of counts for each nucleotide G,A,T,C a numpy
#              array of count values, one per base of the sequence
BaseCounts = namedtuple('BaseCounts', ["G", "A", "T", "C"])


class Cluster:
    """
    The cluster object represents initially, a single molecule bound to a flowcell which is ultimately amplified
    and read from either or both directions.  As such it contains base counts for each position in the original
    molecule's sequence, the sequence(s) called as a result of reading those base counts in either or both directions,
    and the sample barcodes (both 5' and 3' ends).  A collection of these objects forms a cluster_packet.
    """
    # TODO Some coding will change it/when indels are introduced.  They won't easily drop in here.

    next_cluster_id = 1  # Static variable for creating increasing cluster id's

    MIN_ASCII = 33
    MAX_QUALITY = 41

    def __init__(self, run_id, cluster_id, molecule, lane, coordinates, *, molecule_count=1, diameter=0,
                 called_sequences=None, called_barcode=None, quality_scores=None,
                 read_starts = None, read_cigars = None, base_counts=None,
                 forward_is_5_prime=True):
        """
        The constructor contains many attributes, the values of which, may be unknown at the time of instantiation.
        But the object is serializable via custom methods and a serialized version will often contain values for
        these attributes.  So they are in the parameter list here (with defaults) to simplify deserialization.
        :param run_id: The id of the beers run - used in the flowcell header
        :param cluster_id: The unique id of this cluster
        :param molecule: The original molecule object from which the cluster is derived.
        :param lane: The flowcell lane where this cluster is found
        :param coordinates: The other flowcell coordinates (tile, x, y) identifying the cluster's location in the
         flowcell
        :param molecule_count:  The total number of molecules contained within the cluster (starts a 1 and
         geometrically increases with each bridge amplification cycle
        :param diameter: The physical space take up by the cluster (currently unused)
        :param called_sequences: An array of sequences that result from single or paired end reads.  The first
        sequence in the array is the forward read.  The second is the reverse read.  Which directions forward and
        reverse refer to are determined by the 'forward_is_5_prime' attribute.
        :param called_barcode: A string representing the read 5' barcode + 3' barcode - used in the flowcell header
        :param quality_scores: An array of quality score sequences corresponding to the called sequences array.
        :param read_starts: An array of alignment start position, one for each read, aligning them to the reference genome
        :param read_cigars: An array of alignment cigar strings, one for each read, aligning them to the reference genome
        :param base_counts: An array of counts for each of the 4 nt bases for each position in the original molecule.
        The initial base counts array exactly matches the original molecule sequence
        :param forward_is_5_prime: The 5' direction is forward if true.  Otherwise, the 3' direction is forward.
        """
        self.lane = lane
        self.coordinates = coordinates
        self.cluster_id = cluster_id
        self.run_id = run_id
        self.molecule = molecule
        self.diameter = diameter
        self.molecule_count = molecule_count
        self.forward_is_5_prime = forward_is_5_prime
        self.quality_scores = quality_scores or []
        self.called_sequences = called_sequences or []
        self.called_barcode = called_barcode
        self.read_starts = read_starts or []
        self.read_cigars = read_cigars or []
        if base_counts:
            self.base_counts = base_counts
        else:
            # From sequence, we start with 1 count for each base in the sequence
            encoded = np.frombuffer(molecule.sequence.encode("ascii"), dtype='uint8')
            counts = {nt: (encoded == ord(nt)).astype('int32') for nt in "ACGT"}
            self.base_counts = BaseCounts(counts['G'], counts['A'], counts['T'], counts['C'])
        self.forward_is_5_prime = forward_is_5_prime

    def assign_coordinates(self, coordinates):
        """
        Assign a lane coordinates object to the cluster.
        :param coordinates:  a tuple (tile, x, y)
        """
        self.coordinates = coordinates

    def generate_fasta_header(self, direction):
        """
        Assembles the cluster data into a FASTQ header for the given direction (forward or reverse)
        :param direction: The direction for which the header is created (using convention as described in the
        constants module).
        """
        if self.forward_is_5_prime and direction == CONSTANTS.DIRECTION_CONVENTION[0]:
            self.header = f"@{self.encode_sequence_identifier()}\t1:N:0:{self.called_barcode}"
        else:
            self.header = f"@{self.encode_sequence_identifier()}\t2:N:0:{self.called_barcode}"

    def encode_sequence_identifier(self):
        """
        Creates a sequence identifier for the cluster based upon the information contained in the cluster.  The
        instrument is BEERS.  The run id is the identifier set by the user when running the software.  The lane and the
        coordinates define the location on the flowcell.  The flowcell id is just filler for now.
        :return: The sequence identifier as defined for the FASTQ format used here.
        """
        # TODO the 1 is a placeholder for flowcell.  What should we do with this?
        tile, x, y = self.coordinates
        return f"BEERS:{self.run_id}:1:{self.lane}:{tile}:{x}:{y}:{self.molecule.molecule_id}"

    def set_forward_direction(self, forward_is_5_prime):
        """
        Sets which direction (5' or 3') is defined as the forward direction.
        :param forward_is_5_prime: True if the 5' direction is forward and false otherwise.
        """
        self.forward_is_5_prime = forward_is_5_prime

    def read(self, read_length, paired_ends, i5_start, i5_length, i7_start, i7_length, read_start_5_prime, read_start_3_prime):
        """
        Called by the sequencing step in the sequence pipeline to read the cluster in one or both directions.  The
        barcodes on both ends are read as well here.
        :param read_length: The number of bases to be read from either end
        :param paired_ends: Whether only the forward direction is read or both forward and reverse directions are read.
        :param i5_start: 1-based index of first base of i5 barcode from 5' end
        :param i5_length: length of the i5 barcode
        :param i7_start: 1-based index of first base of i7 barcode from 3' end
        :param i7_length: length of the i5 barcode
        :param read_start_5_prime: 1-based index of first base to read (from 5' end)
        :param read_start_3_prime: 1-based index of first base to read (from 3' end)
        reading.
        """
        self.read_barcode(i5_start, i5_length, i7_start, i7_length)
        if self.forward_is_5_prime:
            self.read_in_5_prime_direction(read_length, read_start_5_prime)
            if paired_ends:
                self.read_in_3_prime_direction(read_length, read_start_3_prime)
        else:
            self.read_in_3_prime_direction(read_length, read_start_3_prime)
            if paired_ends:
                self.read_in_5_prime_direction(read_length, read_start_5_prime)

    def read_over_range(self, range_start, range_end):
        """
        Accumulates the called bases and the associated quality scores over the range given (going 5' to 3').  The
        called base at any position is the most numerous base at that position.  In the event of a tie, an 'N' is
        called instead.  The quality score is calculated as a Phred quality score and encoded using the Sanger format.
        :param range_start:  The starting position (from the 5' end), 0-based
        :param range_end: The ending position (exlusive, from the 5' end), 0-based
        """
        with closing(StringIO()) as called_bases:
            with closing(StringIO()) as quality_scores:
                for position in range(range_start, range_end):
                    base_counts = list(zip('GATC', self.get_base_counts_by_position(position)))
                    max_base_count = max(base_counts, key=lambda base_count: base_count[1])
                    number_max_values = len([base_count[0] for base_count in base_counts
                                            if base_count[1] == max_base_count[1]])
                    if number_max_values > 1:
                        called_bases.write('N')
                        quality_scores.write(str(chr(Cluster.MIN_ASCII)))
                    else:
                        other_bases_total_count = sum(count for base, count in base_counts if base != max_base_count[0])
                        prob = other_bases_total_count / self.molecule_count
                        quality_value = min(Cluster.MAX_QUALITY, math.floor(-10 * math.log10(prob))) \
                                        if prob != 0 else Cluster.MAX_QUALITY
                        called_bases.write(max_base_count[0])
                        quality_scores.write(str(chr(quality_value + Cluster.MIN_ASCII)))
                read_start, read_cigar, new_strand = beers_utils.cigar.chain(
                    range_start + 1, # range_start is 0-based, alignments are always 1 based
                    f"{range_end - range_start}M",
                    "+",#TODO: do both 5' and 3' reads use + strand?
                    self.molecule.source_start,
                    self.molecule.source_cigar,
                    self.molecule.source_strand,
                )
                bases = called_bases.getvalue()
                qual = quality_scores.getvalue()
                return bases, qual, read_start, read_cigar

    def read_barcode(self, i5_start, i5_length, i7_start, i7_length):
        """
        Read the barcode information for each end of the cluster using the starting position and length of each barcode
        provided.  The barcode is saved as the 5' portion followed the reverse complemented 3' portion, delimited by
        a plus sign.
        :param i5_start: 1-based index of first base of i5 barcode from 5' read
        :param i5_length: length of the i5 barcode
        :param i7_start: 1-based index of first base of i7 barcode from 3' read
        :param i7_length: length of the i7 barcode
        """
        barcode_5, _, _, _ = self.read_over_range(i5_start - 1, i5_start + i5_length - 1)
        i7_start = len(self.molecule.sequence) - (i7_start - 1) - i7_length
        barcode_3,  _, _, _ = self.read_over_range(i7_start, i7_start + i7_length)
        self.called_barcode = f"{barcode_5}+{GeneralUtils.create_complement_strand(barcode_3)}"

    def read_in_5_prime_direction(self, read_length, read_start):
        """
        Reads the 5' portion of the molecule (skipping over the adapter sequence) out to the given read length.  The
        resulting bases and quality scores are appended to the lists containing called sequences and quality scores,
        respectively.
        :param read_length:  number of bases in from the molecule 5' end to be read
        :param read_start: how far to skip over the adapters to the read (from 5' end)
        """
        range_start = read_start - 1
        range_end = range_start + read_length
        called_bases, quality_scores, read_start, read_cigar = self.read_over_range(range_start, range_end)
        self.quality_scores.append(quality_scores)
        self.called_sequences.append(called_bases)
        self.read_starts.append(read_start)
        self.read_cigars.append(read_cigar)

    def read_in_3_prime_direction(self, read_length, read_start):
        """
        Reads the 3' portion of the molecule (skipping over the adapter sequence) out to the given read length.  The
        resulting bases and quality scores are appended to the lists containing called sequences and quality scores,
        respectively.  Note that the called bases are reverse complemented and the quality scores reversed prior to
        being appended to the aforementioned lists.
        :param read_length:  number of bases in from the molecule 3' end to be read
        :param read_start: how far to skip over the adapters to the read (from 3' end)
        """
        range_start = max(len(self.molecule.sequence) - read_length - (read_start - 1), 0)
        range_end = range_start + read_length
        called_bases, quality_scores, read_start, read_cigar = self.read_over_range(range_start, range_end)
        quality_scores = quality_scores[::-1]
        called_bases = GeneralUtils.create_complement_strand(called_bases)
        self.quality_scores.append(quality_scores)
        self.called_sequences.append(called_bases)
        self.read_starts.append(read_start)
        self.read_cigars.append(read_cigar)

    def get_base_counts_by_position(self, index):
        """
        Conveninece method to gather base counts by index (position)
        :param index: index (position) of interest
        :return: tuple of base counts in GATC order.
        """
        return (self.base_counts.G[index],
                self.base_counts.A[index],
                self.base_counts.T[index],
                self.base_counts.C[index])

    def __str__(self):
        """
        String representation of a cluster that may be displayed when a cluster is printed.  This may not be a
        complete representation.  Depends on what is useful for debugging.
        :return: string representation.
        """
        header = f"run id: {self.run_id}, cluster_id: {self.cluster_id}, molecule_id: {self.molecule.molecule_id}, " \
                 f"molecule_count: {self.molecule_count}, lane: {self.lane}, coordinates: {self.coordinates}\n"
        for index in range(len(self.called_sequences)):
            header += f"called sequence: {self.called_sequences[index]}\n"
            header += f"quality score: {self.quality_scores[index]}\n"
        with closing(StringIO()) as output:
            output.write("pos\tG\tA\tT\tC\torig\n")
            for index in range(len(self.molecule.sequence)):
                output.write(f"{index}\t")
                [output.write(f"{base_count}\t") for base_count in self.get_base_counts_by_position(index)]
                output.write(f"{self.molecule.sequence[index]}\n")
            return header + output.getvalue()

    def serialize(self):
        """
        Method to serialize data into a string.  The initial line contains single value attrbitutes and is preceded
        by a '#'.  The second line contains serialized versions of the coordinates and the molecules and is also
        preceded by a '#'.  Next, may be 0 to 2 lines preceded by '##', each containing a called sequence and
        a called quality score.  The remaining line have no preceding hashes but list the base counts for each
        position on the molecule.
        :return: The serialized string output.
        """
        tile, x, y = self.coordinates
        output = f"#{self.cluster_id}\t{self.run_id}\t{self.molecule_count}\t{self.diameter}\t{self.lane}\t" \
                 f"{self.forward_is_5_prime}\t{self.called_barcode}\n"
        output += f"#{tile}\t{x}\t{y}\n#{self.molecule.serialize()}\n"
        for index in range(len(self.called_sequences)):
            output += f"##{self.called_sequences[index]}\t{self.quality_scores[index]}\t{self.read_starts[index]}\t{self.read_cigars[index]}\n"
        with closing(StringIO()) as counts:
            for index in range(len(self.molecule.sequence)):
                counts.write("\t".join([str(base_count) for base_count in self.get_base_counts_by_position(index)]))
                counts.write("\n")
            output += counts.getvalue()
        return output

    @staticmethod
    def deserialize(data):
        """
        Method to re-render the result of the serialize method above back into a complete cluster object (without
        losing or scrambling state)
        :param data: The string data containing the serialized version of the object
        :return The Cluster object created from the data provided.
        """
        data = data.rstrip('\n')
        called_sequences = []
        quality_scores = []
        read_starts = []
        read_cigars = []
        g_counts = []
        a_counts = []
        t_counts = []
        c_counts = []
        for line_number, line in enumerate(data.split("\n")):
            if line.startswith("##"):
                called_sequence, quality_score, read_start, read_cigar = line[2:].rstrip('\n').split("\t")
                called_sequences.append(called_sequence)
                quality_scores.append(quality_score)
                read_starts.append(int(read_start))
                read_cigars.append(read_cigar)
            elif line.startswith("#"):
                if line_number == 0:
                    cluster_id, run_id, molecule_count, diameter, lane, forward_is_5_prime, called_barcode \
                        = line[1:].rstrip('\n').split("\t")
                if line_number == 1:
                    coordinates = tuple(int(x) for x in (line[1:].rstrip('\n').split("\t")))
                if line_number == 2:
                    molecule = Molecule.deserialize(line[1:].rstrip('\n'))
            else:
                g_count, a_count, t_count, c_count = line.rstrip('\n').split("\t")
                g_counts.append(int(g_count))
                a_counts.append(int(a_count))
                t_counts.append(int(t_count))
                c_counts.append(int(c_count))
        base_counts = BaseCounts(g_counts, a_counts, t_counts, c_counts)
        return Cluster(
                run_id = int(run_id),
                cluster_id = cluster_id,
                molecule = molecule,
                lane = int(lane),
                coordinates = coordinates,
                molecule_count = int(molecule_count),
                diameter = int(diameter),
                called_sequences = called_sequences,
                called_barcode = called_barcode,
                quality_scores = quality_scores,
                base_counts = base_counts,
                read_starts = read_starts,
                read_cigars = read_cigars,
                forward_is_5_prime = forward_is_5_prime
            )
