from collections import namedtuple
import functools
import numpy as np
import scipy.special
import scipy.stats
from io import StringIO
import math
from contextlib import closing
from beers_utils.general_utils import GeneralUtils
from beers_utils.molecule import Molecule
from beers_utils.constants import CONSTANTS
import beers_utils.cigar


BASE_ORDER = ["A", "C", "G", "T"]


class Cluster:
    """
    The cluster object represents initially, a single molecule bound to a flowcell which is ultimately amplified
    and read from either or both directions.  As such it contains base counts for each position in the original
    molecule's sequence, the sequence(s) called as a result of reading those base counts in either or both directions,
    and the sample barcodes (both 5' and 3' ends).  A collection of these objects forms a cluster_packet.
    """
    # TODO Some coding will change it/when indels are introduced.  They won't easily drop in here.

    next_cluster_id = 1  # Static variable for creating increasing cluster id's

    # Sequence-by-synthesis parameters
    PHRED_MIN_ASCII = 33
    PHRED_MAX_QUALITY = 41
    MAX_SKIPS = 10
    MAX_DROPS = 10
    EPSILON = 100 # noise standard deviation
    # Cross-talk between flourescence channels
    # chosen arbitrarily TODO: find a realistic table
    CROSS_TALK = np.array(
        [[1.0, 0.3, 0.2, 0.0],
         [0.1, 1.0, 0.2, 0.1],
         [0.2, 0.1, 1.0, 0.3],
         [0.3, 0.1, 0.1, 1.0]]
        )

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
        if base_counts is not None:
            self.base_counts = base_counts
        else:
            # From sequence, we start with 1 count for each base in the sequence
            encoded = np.frombuffer(molecule.sequence.encode("ascii"), dtype='uint8')
            self.base_counts = np.array([encoded == ord(nt) for nt in BASE_ORDER]).astype(int)
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

    def read_over_range(self, range_start, range_end, skip_rate = 0 , drop_rate = 0):
        """
        Accumulates the called bases and the associated quality scores over the range given (going 5' to 3').  The
        called base at any position is the most numerous base at that position.  In the event of a tie, an 'N' is
        called instead.  The quality score is calculated as a Phred quality score and encoded using the Sanger format.
        :param range_start:  The starting position (from the 5' end), 0-based
        :param range_end: The ending position (exlusive, from the 5' end), 0-based
        :param skip_rate: rate per base at which two bases are added instead of one in the sequence-by-synthesis
        :param drop_rate: rate per base at which no base is added instead of one in the sequence-by-synthesis
        """
        base_array = np.array([ord(x) for x in BASE_ORDER], dtype="uint8")
        with closing(StringIO()) as called_bases:
            with closing(StringIO()) as quality_scores:

                base_counts = self.base_counts[:, range_start:range_end]
                padded_base_counts = np.concatenate(
                        (   np.zeros((4,self.MAX_DROPS)),
                            self.base_counts,
                            np.zeros((4,self.MAX_SKIPS))
                        ),
                        axis=1)
                read_len = range_end - range_start

                #TODO: this all depends upon the direction of the read, forward or reverse?
                #TODO: this should all be done at the level of a lane or tile or something
                #      various parameters need to be estimated across all clusters

                skip_rate = 0.003
                drop_rate = 0.003

                # Add skips (aka prephasing) where bases are added too quickly in some strands
                # This causes bases to be read too early
                # We assume that each molecule skips at most MAX_SKIPS times
                # And we approximate skips and drops without simulating every molecule
                # individually by just using average molecule counts

                # Each base is equally likely to give a skip
                skips = np.random.random(size = (self.molecule_count, range_end - range_start)) < skip_rate
                num_skips_so_far = np.minimum(np.cumsum(skips, axis=1), self.MAX_SKIPS)
                # Count number of molecules that have had exactly k skips at the nth base, for k= 0,1,...MAX_SKIPS
                num_mols_skipped = (num_skips_so_far[:,:,None] == np.arange(0, self.MAX_SKIPS+1)[None,None,:]).sum(axis=0)
                frac_skipped = num_mols_skipped / self.molecule_count

                # Each base is equally likely to give a drop
                drops = np.random.random(size = (self.molecule_count, range_end - range_start)) < drop_rate
                num_drops_so_far = np.minimum(np.cumsum(drops, axis=1), self.MAX_SKIPS)
                # Count number of molecules that have had exactly k drops at the nth base, for k= 0,1,...MAX_DROPS
                num_mols_dropped = (num_drops_so_far[:,:,None] == np.arange(0, self.MAX_DROPS+1)[None,None,:]).sum(axis=0)
                frac_dropped = num_mols_dropped / self.molecule_count

                # 'smear' base counts according to the number of skips
                smeared_base_counts = sum( padded_base_counts[:, self.MAX_DROPS + range_start + skip : self.MAX_DROPS + range_end + skip] * frac_skipped[:,skip]
                                                for skip in range(1, self.MAX_SKIPS+1))
                # ...and by the number of drops
                smeared_base_counts = sum( padded_base_counts[:, self.MAX_DROPS + range_start - drop : self.MAX_DROPS + range_end - drop] * frac_dropped[:,drop]
                                                for drop in range(1, self.MAX_DROPS+1))
                # and the ones that neither drop nor skip
                smeared_base_counts = padded_base_counts[:, range_start:range_end] * (1 -  (1 - frac_skipped[:, 0]) - (1 - frac_dropped[:, 0]))
                
                # Flourescence comes from these after including cross talk between the different colors
                flourescence = self.CROSS_TALK @ (smeared_base_counts + self.EPSILON * np.random.normal(size=smeared_base_counts.shape))

                # Find the brightest base at each position
                #brightest_base = np.argmax(flourescence, axis=0) # bases as nums 0,1,2,3
                #read_seq = base_array[brightest_base].tobytes().decode() # as string "ACGT..."

                # Approximate quality scores from the ratio of maximum base count
                # to the second most common base
                #ordered_flourescence = np.sort(flourescence, axis=0)
                #highest_flourescence = ordered_flourescence[-1,:]
                #second_highest_flourescence = ordered_flourescence[-2,:]
                ##TODO: these probabilites are super high
                #prob = second_highest_flourescence / (second_highest_flourescence + highest_flourescence)
                ##combos_highest = scipy.special.comb(highest_flourescence + second_highest_flourescence, highest_flourescence)
                ##combos_second = scipy.special.comb(highest_flourescence + second_highest_flourescence, second_highest_flourescence)
                ##prob2 = scipy.stats.binom(self.molecule_count, p=0.5).sf(highest_flourescence)
                #prob = np.maximum(prob, np.power(10, -Cluster.PHRED_MAX_QUALITY / 10)) # Force minimum phred score

                ## Approximate the Bustard algorithm for calling bases and quality scores
                # make an 'estimated' cross-talk matrix TODO: should be estimated from an entire tile
                cross_talk_est = self.CROSS_TALK + np.random.normal(size=self.CROSS_TALK.shape)* 0.01
                cross_talk_inv_est = np.linalg.inv(cross_talk_est)
                #TODO: bustard also accounts for phasing/prephasing
                # gives probability of a template terminating at position j after t cycles
                phasing_matrix_inv = get_inv_phasing_matrix(read_len, skip_rate, drop_rate)
                base_counts_est = cross_talk_inv_est @ flourescence @ phasing_matrix_inv
                base_orders = np.argsort(base_counts_est, axis=0)
                called_base = base_orders[-1,:] # bases as nums 0,1,2,3 = ACGT
                second_best_base = base_orders[-2,:] # bases as nums 0,1,2,3 = ACGT
                read_seq = base_array[called_base].tobytes().decode() # as string "ACGT..."
                #TODO: epsilon should be estimated from data and is cycle-dependent
                M = (cross_talk_inv_est * self.EPSILON) * (cross_talk_inv_est * self.EPSILON).T
                highest_base_count = base_counts_est[called_base, np.arange(read_len)]
                second_base_count = base_counts_est[second_best_base, np.arange(read_len)]
                diff = highest_base_count - second_base_count
                contrasts = np.zeros(shape=(read_len, 4, 1))
                contrasts[np.arange(read_len), called_base] = 1
                contrasts[np.arange(read_len), second_best_base] = -1
                sigma = np.sqrt(contrasts.transpose(0,2,1) @ M @ contrasts)
                prob = scipy.stats.norm(0, sigma.flatten()).sf(np.abs(diff)) * 2 # two-tailed

                score = (-10 * np.log10(prob)).astype("uint8")
                score = np.minimum(Cluster.PHRED_MAX_QUALITY, score)
                qual = (Cluster.PHRED_MIN_ASCII + score).tobytes().decode()
                seq = np.array(BASE_ORDER)[called_base].tobytes().decode()

                #print("Flourescence")
                #print(flourescence[called_base, np.arange(flourescence.shape[1])])
                #print(flourescence[second_best_base, np.arange(flourescence.shape[1])])
                #print(prob[:10])
                #print(qual[:10])

                read_start, read_cigar, new_strand = beers_utils.cigar.chain(
                    range_start + 1, # range_start is 0-based, alignments are always 1 based
                    f"{range_end - range_start}M",
                    "+",#TODO: do both 5' and 3' reads use + strand?
                    self.molecule.source_start,
                    self.molecule.source_cigar,
                    self.molecule.source_strand,
                )
                return read_seq, qual, read_start, read_cigar

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
        :return: tuple of base counts in ACGT order.
        """
        return tuple(self.base_counts[:,index])

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
            output.write("pos\tA\tC\tG\tT\torig\n")
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
                a_count, c_count, g_count, t_count = line.rstrip('\n').split("\t")
                a_counts.append(int(a_count))
                c_counts.append(int(c_count))
                g_counts.append(int(g_count))
                t_counts.append(int(t_count))
        base_counts = np.array((a_counts, c_counts, g_counts, t_counts))
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

def hypergeometric(counts, n_samples):
    '''
    Vectorized version of np.random.hypergeometric

    Draw n_sample balls from an urn with counts.shape[0] colors of balls

    :param counts: n-dim array with first dimension giving the counts of 'balls' available of each of counts.shape[0] colors
    :params n_samples: n-1-dim array with number of samples to draw from the urn
    '''

    def hypergeometric_(good_counts, bad_counts, n_samples):
        # same as np.random.hypergeometric but allows n_samples = 0
        out = np.zeros_like(good_counts)
        nonzero = n_samples > 0
        out[nonzero] = np.random.hypergeometric(good_counts[nonzero], bad_counts[nonzero], n_samples[nonzero])
        return out

    assert counts.shape[1:] == n_samples.shape

    out = np.empty_like(counts)

    total_counts = counts.sum(axis=0)
    assert (n_samples <= total_counts).all()

    # Iterate through the 'colors', drawing some balls of each color and passing
    # the remaining draws to the remaining colors
    n_samples = n_samples.copy()
    for i in range(counts.shape[0]-1):
        color_counts = counts[i]
        draws_from_this_color = hypergeometric_(color_counts, total_counts - color_counts, n_samples)
        out[i] = draws_from_this_color
        n_samples -= draws_from_this_color
        total_counts -= color_counts
    out[-1] = n_samples # All remaining go to last color
    return out

@functools.lru_cache(maxsize=None)
def get_inv_phasing_matrix(read_len, skip_rate, drop_rate):
    ''' Computes the inverse of Q, the phasing matrix

    The (j,t) entry of Q gives the probability of a template 
    terminating at position j after t cycles

    Cached for speed-ups. Since reads are all the same length, they all get the same
    matrix and caching is very useful.
    '''
    phasing_matrix = np.array([[(1 -  skip_rate - drop_rate) if j == t else 
                                    (skip_rate**(j - t)*(1-skip_rate) if j > t else
                                     drop_rate**(t - j)*(1-drop_rate))
                                for t in range(read_len)]
                                for j in range(read_len)])
    # Strangely, computing this inverse is very slow on cluster (~3 seconds)
    # even for small size (100x100), though very fast on my laptop (~3ms)
    phasing_mat_inv = np.linalg.inv(phasing_matrix)
    return phasing_mat_inv
