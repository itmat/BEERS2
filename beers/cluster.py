from collections import namedtuple
import numpy as np
from io import StringIO
import math
from contextlib import closing
from beers.utilities.general_utils import GeneralUtils
from beers.flowcell_lane import LaneCoordinates
from beers_utils.molecule import Molecule
from beers_utils.constants import CONSTANTS


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

    def __init__(self, run_id, cluster_id, molecule, lane, coordinates, molecule_count=1, diameter=0,
                 called_sequences=None, called_barcode=None, quality_scores=None, base_counts=None,
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
        encoded = np.frombuffer(molecule.sequence.encode("ascii"), dtype='uint8')
        counts = {nt: (encoded == ord(nt)).astype('int32') for nt in "ACGT"}
        if base_counts:
            self.base_counts = base_counts
        else:
            self.base_counts = BaseCounts(counts['G'], counts['A'], counts['T'], counts['C'])
        self.forward_is_5_prime = forward_is_5_prime

    def assign_coordinates(self, coordinates):
        """
        Assign a lane coordinates object to the cluster.
        :param coordinates: An object of the LaneCoordinates class
        """
        self.coordinates = coordinates

    def generate_fasta_header(self, direction):
        """
        Assembles the cluster data into a FASTQ header for the given direction (forward or reverse)
        :param direction: The direction for which the header is created (using convention as described in the
        constants module).
        """
        if self.forward_is_5_prime and direction == CONSTANTS.DIRECTION_CONVENTION[0]:
            self.header = f"{self.encode_sequence_identifier()}\t1:N:0:{self.called_barcode}"
        else:
            self.header = f"{self.encode_sequence_identifier()}\t2:N:0:{self.called_barcode}"

    def encode_sequence_identifier(self):
        """
        Creates a sequence identifier for the cluster based upon the information contained in the cluster.  The
        instrument is BEERS.  The run id is the identifier set by the user when running the software.  The lane and the
        coordinates define the location on the flowcell.  The flowcell id is just filler for now.
        :return: The sequence identifier as defined for the FASTQ format used here.
        """
        # TODO the 1 is a placeholder for flowcell.  What should we do with this?
        return f"@BEERS:{self.run_id}:1:{self.lane}:{self.coordinates.tile}:{self.coordinates.x}:{self.coordinates.y}"

    def set_forward_direction(self, forward_is_5_prime):
        """
        Sets which direction (5' or 3') is defined as the forward direction.
        :param forward_is_5_prime: True if the 5' direction is forward and false otherwise.
        """
        self.forward_is_5_prime = forward_is_5_prime

    def read(self, read_length, paired_ends, barcode_data, adapter_sequences):
        """
        Called by the sequencing step in the sequence pipeline to read the cluster in one or both directions.  The
        barcodes on both ends are read as well here.
        :param read_length: The number of bases to be read from either end
        :param paired_ends: Whether only the forward direction is read or both forward and reverse directions are read.
        :param barcode_data: The start position and length of each barcode so that they may be located.
        :param adapter_sequences: The sequences of the adapters on each end.  Used to skip over the adapters when
        reading.
        """
        self.read_barcode(barcode_data)
        if self.forward_is_5_prime:
            self.read_in_5_prime_direction(read_length, adapter_sequences[0])
            if paired_ends:
                self.read_in_3_prime_direction(read_length, adapter_sequences[1])
        else:
            self.read_in_3_prime_direction(read_length, adapter_sequences[1])
            if paired_ends:
                self.read_in_5_prime_direction(read_length, adapter_sequences[0])

    def read_over_range(self, range_start, range_end):
        """
        Accumulates the called bases and the associated quality scores over the range given (going 5' to 3').  The
        called base at any position is the most numerous base at that position.  In the event of a tie, an 'N' is
        called instead.  The quality score is calculated as a Phred quality score and encoded using the Sanger format.
        :param range_start:  The starting position (from the 5' end)
        :param range_end: The ending position (exlusive, from the 5' end)
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
                return called_bases.getvalue(), quality_scores.getvalue()

    def read_barcode(self, barcode_data):
        """
        Read the barcode information for each end of the cluster using the starting position and length of each barcode
        provided.  The barcode is saved as the 5' portion followed the reverse complemented 3' portion, delimited by
        a plus sign.
        :param barcode_data: An array containing the starting positions and lengths of each barcode.
        """
        range_start = barcode_data[0]
        range_end = barcode_data[0] + barcode_data[1]
        barcode_5, _ = self.read_over_range(range_start, range_end)
        range_end = len(self.molecule.sequence) - barcode_data[2]
        range_start = range_end - barcode_data[3]
        barcode_3,  _ = self.read_over_range(range_start, range_end)
        self.called_barcode = f"{barcode_5}+{GeneralUtils.create_complement_strand(barcode_3)}"

    def read_in_5_prime_direction(self, read_length, adapter_sequence):
        """
        Reads the 5' portion of the molecule (skipping over the adapter sequence) out to the given read length.  The
        resulting bases and quality scores are appended to the lists containing called sequences and quality scores,
        respectively.
        :param read_length:  number of bases in from the molecule 5' end to be read
        :param adapter_sequence: the 5' adapter sequence to skip over.
        """
        range_start = len(adapter_sequence)
        range_end = range_start + read_length
        called_bases, quality_scores = self.read_over_range(range_start, range_end)
        self.quality_scores.append(quality_scores)
        self.called_sequences.append(called_bases)

    def read_in_3_prime_direction(self, read_length, adapter_sequence):
        """
        Reads the 3' portion of the molecule (skipping over the adapter sequence) out to the given read length.  The
        resulting bases and quality scores are appended to the lists containing called sequences and quality scores,
        respectively.  Note that the called bases are reverse complemented and the quality scores reversed prior to
        being appended to the aforementioned lists.
        :param read_length:  number of bases in from the molecule 3' end to be read
        :param adapter_sequence: the 3' adapter sequence to skip over.
        """
        range_end = len(self.molecule.sequence) - len(adapter_sequence)
        range_start = range_end - read_length
        called_bases, quality_scores = self.read_over_range(range_start, range_end)
        quality_scores = quality_scores[::-1]
        called_bases = GeneralUtils.create_complement_strand(called_bases)
        self.quality_scores.append(quality_scores)
        self.called_sequences.append(called_bases)

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
        output = f"#{self.cluster_id}\t{self.run_id}\t{self.molecule_count}\t{self.diameter}\t{self.lane}\t" \
                 f"{self.forward_is_5_prime}\t{self.called_barcode}\n"
        output += f"#{self.coordinates.serialize()}\n#{self.molecule.serialize()}\n"
        for index in range(len(self.called_sequences)):
            output += f"##{self.called_sequences[index]}\t{self.quality_scores[index]}\n"
        with closing(StringIO()) as counts:
            for index in range(len(self.molecule.sequence)):
                [counts.write(f"{base_count}\t") for base_count in self.get_base_counts_by_position(index)]
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
        data = data.rstrip()
        called_sequences = []
        quality_scores = []
        g_counts = []
        a_counts = []
        t_counts = []
        c_counts = []
        for line_number, line in enumerate(data.split("\n")):
            if line.startswith("##"):
                called_sequence, quality_score = line[2:].rstrip().split("\t")
                called_sequences.append(called_sequence)
                quality_scores.append(quality_score)
            elif line.startswith("#"):
                if line_number == 0:
                    cluster_id, run_id, molecule_count, diameter, lane, forward_is_5_prime, called_barcode \
                        = line[1:].rstrip().split("\t")
                if line_number == 1:
                    coordinates = LaneCoordinates.deserialize(line[1:].rstrip())
                if line_number == 2:
                    molecule = Molecule.deserialize(line[1:].rstrip())
            else:
                g_count, a_count, t_count, c_count = line.rstrip().split("\t")
                g_counts.append(int(g_count))
                a_counts.append(int(a_count))
                t_counts.append(int(t_count))
                c_counts.append(int(c_count))
        base_counts = BaseCounts(g_counts, a_counts, t_counts, c_counts)
        return Cluster(int(run_id), cluster_id, molecule, int(lane), coordinates, int(molecule_count), int(diameter),
                       called_sequences, called_barcode, quality_scores, base_counts, forward_is_5_prime)
