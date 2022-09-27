from io import StringIO
import math
from typing import Optional, ClassVar
from dataclasses import dataclass, field
from contextlib import closing
import numpy as np
from beers_utils.general_utils import GeneralUtils
from beers_utils.molecule import Molecule
from beers_utils.constants import CONSTANTS
import beers_utils.cigar


BASE_ORDER = ["A", "C", "G", "T"]



@dataclass
class Cluster:
    """
    The cluster object represents initially, a single molecule bound to a flowcell which is ultimately amplified
    and read from either or both directions.  As such it contains base counts for each position in the original
    molecule's sequence, the sequence(s) called as a result of reading those base counts in either or both directions,
    and the sample barcodes (both 5' and 3' ends).  A collection of these objects forms a cluster_packet.

    Attributes
    ----------
    cluster_id:
        The unique id of this cluster
    molecule:
        The original molecule object from which the cluster is derived.
    lane:
        The flowcell lane where this cluster is found
    coordinates:
        The other flowcell coordinates (tile, x, y) identifying the cluster's location in the flowcell
    molecule_count:
        The total number of molecules contained within the cluster (starts a 1 and geometrically increases
        with each bridge amplification cycle
    diameter:
        The physical space take up by the cluster (currently unused)
    called_sequences:
        An array of sequences that result from single or paired end reads.  The first
        sequence in the array is the forward read.  The second is the reverse read.
    called_barcode:
        A string representing the read 5' barcode + 3' barcode - used in the flowcell header
    quality_scores:
        An array of quality score sequences corresponding to the called sequences array.
    read_starts:
        An array of alignment start position, one for each read, aligning them to the reference genome
    read_cigars:
        An array of alignment cigar strings, one for each read, aligning them to the reference genome
    read_strands:
        An array of alignment strand strings, one for each read, aligning them to the reference genome
    base_counts:
        An 4xn array of counts for each of the 4 nt bases for each position in the original molecule.
        In the order of A, C, G, T counts.  The initial base counts array exactly matches the original
        molecule sequence.
    """

    cluster_id: int
    molecule: Molecule
    lane: int
    coordinates: tuple[int, int, int]
    molecule_count: int =1
    diameter: int = 0
    called_sequences: list[str] = field(default_factory=list)
    called_barcode: Optional[str] = None
    quality_scores: list[str] = field(default_factory=list)
    read_starts: list[int] = field(default_factory=list)
    read_cigars: list[str] = field(default_factory=list)
    read_strands: list[str] = field(default_factory=list)
    base_counts: Optional[np.ndarray] = None

    next_cluster_id: ClassVar[float] = 1  # Static variable for creating increasing cluster id's

    def initialize_base_counts(self):
        """
        Generates basecounts from the starting molecule
        """
        # From sequence, we start with 1 count for each base in the sequence
        encoded = np.frombuffer(self.molecule.sequence.encode("ascii"), dtype='uint8')
        self.base_counts = np.array([encoded == ord(nt) for nt in BASE_ORDER]).astype(int)
        self.molecule_count = 1

    def assign_coordinates(self, coordinates: tuple[int, int, int]):
        """
        Assign a lane coordinates object to the cluster.

        Parameters
        ----------
        coordinates:
            (tile, x, y)
        """
        self.coordinates = coordinates

    def generate_fasta_header(self, direction: int):
        """
        Assembles the cluster data into a FASTQ header for the given direction (forward or reverse)

        Parameters
        ---------
        direction: int, either 1 or 2
            Direction of the read. 1 is the 'forward' read, 2 is the 'reverse' read.
            If single ended reads, just use 1.
        """
        if direction == CONSTANTS.DIRECTION_CONVENTION[0]:
            self.header = f"@{self.encode_sequence_identifier()}\t1:N:0:{self.called_barcode}"
        else:
            self.header = f"@{self.encode_sequence_identifier()}\t2:N:0:{self.called_barcode}"

    def encode_sequence_identifier(self) -> str:
        """
        Creates a sequence identifier for the cluster based upon the information contained in the cluster.  The
        instrument is BEERS.  The run id is the identifier set by the user when running the software.  The lane and the
        coordinates define the location on the flowcell.  The flowcell id is just filler for now.

        Returns
        -------
        The sequence identifier as defined for the FASTQ format used here.
        """
        # TODO the 1 is a placeholder for flowcell.  What should we do with this?
        tile, x, y = self.coordinates
        return f"BEERS:1:{self.lane}:{tile}:{x}:{y}:{self.molecule.molecule_id}"


    def get_base_counts_by_position(self, index: int) -> tuple[int, int, int, int]:
        """
        Conveninece method to gather base counts by index (position)

        Parameters
        ---------
        index:
            index (position) of interest

        Returns
        ------
        ndarray
            tuple of base counts in ACGT order.
        """
        return tuple(self.base_counts[:,index])  # type: ignore

    def __str__(self) -> str:
        """
        String representation of a cluster that may be displayed when a cluster is printed.  This may not be a
        complete representation.  Depends on what is useful for debugging.

        Returns string representation.
        """
        header = f"cluster_id: {self.cluster_id}, molecule_id: {self.molecule.molecule_id}, " \
                 f"molecule_count: {self.molecule_count}, lane: {self.lane}, coordinates: {self.coordinates}\n"
        for index in range(len(self.called_sequences)):
            header += f"called sequence: {self.called_sequences[index]}\n"
            header += f"quality score: {self.quality_scores[index]}\n"
        with closing(StringIO()) as output:
            output.write("pos\tA\tC\tG\tT\torig\n")
            for index in range(len(self.molecule.sequence)):
                output.write(f"{index}\t")
                [output.write(f"{base_count}\t") for base_count in self.get_base_counts_by_position(index)]
                output.write(f"{self.molecule.sequence[index]}\n")
            return header + output.getvalue()

    def serialize(self) -> str:
        """
        Method to serialize data into a string.  The initial line contains single value attrbitutes and is preceded
        by a '#'.  The second line contains serialized versions of the coordinates and the molecules and is also
        preceded by a '#'.  Next, may be 0 to 2 lines preceded by '##', each containing a called sequence and
        a called quality score.  The remaining line have no preceding hashes but list the base counts for each
        position on the molecule.

        Returns
        -------
        str
            The serialized string output.
        """
        tile, x, y = self.coordinates
        output = f"#{self.cluster_id}\t{self.molecule_count}\t{self.diameter}\t{self.lane}\t" \
                 f"{self.called_barcode}\n"
        output += f"#{tile}\t{x}\t{y}\n#{self.molecule.serialize()}\n"
        for index in range(len(self.called_sequences)):
            output += f"##{self.called_sequences[index]}\t{self.quality_scores[index]}\t{self.read_starts[index]}\t{self.read_cigars[index]}\t{self.read_strands[index]}\n"
        with closing(StringIO()) as counts:
            if self.molecule_count == 1 or self.base_counts is None:
                # Shortcut: we only have the one molecule which is already saved
                # so we indicate that here and don't write out all the counts
                counts.write("None\n")
            else:
                for index in range(len(self.molecule.sequence)):
                    counts.write("\t".join([str(base_count) for base_count in self.get_base_counts_by_position(index)]))
                    counts.write("\n")
                output += counts.getvalue()
        return output

    @staticmethod
    def deserialize(data: str, skip_base_counts: bool =False):
        """
        Method to re-render the result of the serialize method above back into a complete cluster object (without
        losing or scrambling state)

        Parameters
        ----------
        data:
            The string data containing the serialized version of the object
        skip_base_counts:
            If True, don't load the base_count data (for memory efficiency)

        Returns
        -------
        Cluster
            The Cluster object created from the data provided.
        """
        data = data.rstrip('\n')
        called_sequences = []
        quality_scores = []
        read_starts = []
        read_cigars = []
        read_strands = []
        base_counts_not_supplied = False
        g_counts = []
        a_counts = []
        t_counts = []
        c_counts = []
        for line_number, line in enumerate(data.split("\n")):
            if line.startswith("##"):
                called_sequence, quality_score, read_start, read_cigar, read_strand = line[2:].rstrip('\n').split("\t")
                called_sequences.append(called_sequence)
                quality_scores.append(quality_score)
                read_starts.append(int(read_start))
                read_cigars.append(read_cigar)
                read_strands.append(read_strand)
            elif line.startswith("#"):
                if line_number == 0:
                    cluster_id_, molecule_count_, diameter_, lane_, called_barcode \
                        = line[1:].rstrip('\n').split("\t")
                    cluster_id = int(cluster_id_)
                    molecule_count = int(molecule_count_)
                    diameter = int(diameter_)
                    lane = int(lane_)
                if line_number == 1:
                    coordinates = tuple(int(x) for x in (line[1:].rstrip('\n').split("\t")))
                if line_number == 2:
                    molecule = Molecule.deserialize(line[1:].rstrip('\n'))
            elif line == "None" or skip_base_counts:
                base_counts_not_supplied = True
            else:
                a_count, c_count, g_count, t_count = line.rstrip('\n').split("\t")
                a_counts.append(int(a_count))
                c_counts.append(int(c_count))
                g_counts.append(int(g_count))
                t_counts.append(int(t_count))

        if base_counts_not_supplied:
            base_counts = None
        else:
            base_counts = np.array((a_counts, c_counts, g_counts, t_counts))

        return Cluster(
                cluster_id = cluster_id,
                molecule = molecule,
                lane = int(lane),
                coordinates = coordinates, # type: ignore
                molecule_count = int(molecule_count),
                diameter = int(diameter),
                called_sequences = called_sequences,
                called_barcode = called_barcode,
                quality_scores = quality_scores,
                base_counts = base_counts,
                read_starts = read_starts,
                read_cigars = read_cigars,
                read_strands = read_strands,
            )
