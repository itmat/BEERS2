import functools

import numpy as np
import scipy.stats

import warnings

from beers.sequence.sequence_by_synthesis_helper import get_frac_skipped as get_frac_skipped_cython

import beers_utils
import beers_utils.molecule
from beers_utils.general_utils import GeneralUtils
import beers
from  beers.cluster import Cluster
from beers.cluster_packet import ClusterPacket

# Sequence-by-synthesis parameters
PHRED_MIN_ASCII = 33
PHRED_MAX_QUALITY = 41
BASE_ARRAY = np.array([ord(x) for x in beers.cluster.BASE_ORDER], dtype="uint8")

class SequenceBySynthesisStep:
    """
    Sequence by Synthesis Step

    Given a bridge-amplified cluster packet, generates the reads from them.
    Reads are determined by flourescence as bases are added one-by-one.
    We approximate this, allowing for phasing (where some molecules are
    further ahead or behind others in the sequencing steps).

    Flourescence reads are then turned into base calls. Then the values are
    filled in on the Clusters for the reads, in:
    - called_sequences
    - called_barcodes
    - quality_scores
    - read_starts
    - read_cigars
    - read_strands
    Which together give all the read information necessary to make a fastq
    and the ideal true alignment necessary to report as a SAM file with
    true alignment.

    Config Example::

        parameters:
            # Determines which read is the forward and which is the reverse read
            # NOTE: currently only 'true' is implemented
            forward_is_5_prime: true
            # Whether to sequence both ends or just one
            # NOTE: currently only 'true' is implemented
            paired_ends: true
            # Length of the reads to generate, in bases
            read_length: 100
            # The rate of phasing, either forward (skip) or backwards (drop)
            # per base per molecule
            skip_rate: 0.002
            drop_rate: 0.002
    """

    name = "Sequence By Synthesis Step"

    # Sequencing parameters
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

    def __init__(self, step_log_file_path: str, parameters: dict, global_config: dict):
        """
        Initializes the step with a file path to the step log and a dictionary of parameters.  Missing parameters that
        control non-idealized behavior are replaced with default values that produce idealized step behavior.

        Parameters
        ----------
        step_log_file_path:
            location of step logfile
        parameters:
            json-like object of parameters.  Any required parameters not provided are identified by the
            validate method.
        global_config:
            json-like object general parameters, for the overall BEERS run
        """
        self.log_filepath = step_log_file_path
        self.read_length = parameters['read_length']
        self.forward_is_5_prime = parameters.get('forward_is_5_prime', True)
        if self.forward_is_5_prime != True:
            # TODO: support both fwd and rev reads being the first read
            raise NotImplementedError("Sequencing is only supported with forward_is_5_prime = true")
        self.paired_ends = parameters.get('paired_ends', False)
        if self.paired_ends == False:
            #TODO: support single end reads
            raise NotImplementedError("Sequencing is only supported for paired end data at the moment (paired_end = True)")

        # Sequencing error rate
        self.skip_rate = parameters['skip_rate']
        self.drop_rate = parameters['drop_rate']

        # Determine where the barcodes lie
        self.i5_start = len(global_config['resources']['pre_i5_adapter'])
        self.i7_start = len(global_config['resources']['post_i7_adapter'])
        self.post_i5_length = len(global_config['resources']['post_i5_adapter'])
        # Part read after the i7 is the 'pre' adapter, since reading from 3' end
        # NOTE: + 1 for the additional 'A' that is ligated to 3' end before adapter ligation
        self.post_i7_length = len(global_config['resources']['pre_i7_adapter']) + 1
        i5_barcode_lengths = [len(sample['barcodes']['i5']) for sample in  global_config['samples'].values()]
        assert len(set(i5_barcode_lengths)) == 1
        self.i5_length = i5_barcode_lengths[0]
        i7_barcode_lengths = [len(sample['barcodes']['i7']) for sample in  global_config['samples'].values()]
        assert len(set(i7_barcode_lengths)) == 1
        self.i7_length = i7_barcode_lengths[0]

        # Determine where the sequencing starts
        # Read from the start of the i5 barcode, through the rest of the adatpter and into the read
        # TODO: is this actually how the sequencing is done?
        self.forward_read_start = self.i5_start
        self.forward_read_end = self.forward_read_start + self.i5_length + self.post_i5_length + self.read_length
        self.reverse_read_start = self.i7_start
        self.reverse_read_end = self.reverse_read_start + self.i7_length + self.post_i7_length + self.read_length

        print(f"{SequenceBySynthesisStep.name} instantiated")

    def execute(self, cluster_packet: ClusterPacket, rng: np.random.Generator) -> ClusterPacket:
        """
        Execute the Sequence By Synthess step on one packet
        Parameters
        ----------
        cluster_packet:
            The input cluster packet
        rng:
            Random number generator instance

        Returns
        ------
            The updated cluster packet
        """
        print(f"Starting {self.name}")

        with open(self.log_filepath, "wt") as log:

            # make an 'estimated' cross-talk matrix TODO: should be estimated from an entire tile
            cross_talk_est = self.CROSS_TALK + rng.normal(size=self.CROSS_TALK.shape)* 0.01
            cross_talk_inv_est = np.linalg.inv(cross_talk_est)
            # Assume epsilon (flourescence imaging noise) is perfectly known)
            #TODO: epsilon should be estimated from data and can be cycle number dependent
            epsilon_est = self.EPSILON

            for cluster in cluster_packet.clusters:
                # TODO: use self.forward_is_5_prime to know which direction is 'forward'?
                # Perform sequence-by-synthesis and get 'flourescence images'
                forward_flourescence, fwd_start, fwd_cigar, fwd_strand = self.read_flourescence(
                        cluster,
                        self.forward_read_start,
                        self.forward_read_end,
                        direction = "+",
                        rng = rng,
                )
                forward_bases, forward_quality = self.call_bases(forward_flourescence, epsilon_est, cross_talk_inv_est)

                reverse_flourescence, rev_start, rev_cigar, rev_strand = self.read_flourescence(
                        cluster,
                        self.reverse_read_start,
                        self.reverse_read_end,
                        direction = "-",
                        rng = rng,
                )
                reverse_bases, reverse_quality = self.call_bases(reverse_flourescence, epsilon_est, cross_talk_inv_est)

                #Extract barcodes and read sequences
                forward_barcode = forward_bases[:self.i5_length]
                forward_read = forward_bases[self.i5_length + self.post_i5_length:]
                forward_quality = forward_quality[self.i5_length + self.post_i5_length:]

                reverse_barcode = GeneralUtils.create_complement_strand(reverse_bases[:self.i7_length])
                reverse_read = reverse_bases[self.i7_length + self.post_i7_length:]
                reverse_quality = reverse_quality[self.i7_length + self.post_i7_length:]

                cluster.called_sequences = [forward_read, reverse_read]
                cluster.called_barcode = f"{forward_barcode}+{reverse_barcode}"
                cluster.quality_scores = [forward_quality, reverse_quality]

                # Get alignment of read sequences, not including barcodes
                fwd_start, fwd_cigar, fwd_strand = beers_utils.cigar.chain(
                    self.i5_length + self.post_i5_length + 1, # 1-based starts
                    f"{self.read_length}M",
                    "+",
                    fwd_start, fwd_cigar, fwd_strand
                )
                rev_start, rev_cigar, rev_strand = beers_utils.cigar.chain(
                    self.i7_length + self.post_i7_length, # 1-based starts
                    f"{self.read_length}M",
                    "+",
                    rev_start, rev_cigar, rev_strand
                )
                cluster.read_starts = [fwd_start, rev_start]
                cluster.read_cigars = [fwd_cigar, rev_cigar]
                cluster.read_strands = [fwd_strand, rev_strand]

        cluster_packet.clusters = sorted(cluster_packet.clusters, key=lambda cluster: cluster.coordinates)
        return cluster_packet

    @staticmethod
    def validate(parameters: dict, global_config: dict):
        errors = []

        if 'read_length' not in parameters:
            errors.append("Must specify 'read_length'")
        elif not isinstance(parameters['read_length'], int):
            errors.append("'read_length' must be a positive integer")
        elif not parameters['read_length'] > 0:
            errors.append("'read_length' must be a positive")

        if not isinstance(parameters.get('forward_is_5_prime', True), bool):
            errors.append("forward_is_5_prime must be either true (default) or false")

        if not isinstance(parameters.get('paired_ends', True), bool):
            errors.append("paired_ends must be either true (default) or false")

        for var in ['skip_rate', 'drop_rate']:
            if var not in parameters:
                errors.append(f"Must specify {var}")
            elif not isinstance(parameters[var], (float, int)):
                errors.append(f"{var} must be a number")
            elif parameters[var] < 0:
                errors.append(f"{var} must be non-negative")

        # Validate adapters
        resources =  global_config['resources']
        for var in ['pre_i5_adapter', 'post_i7_adapter', 'post_i5_adapter', 'pre_i5_adapter']:
            if var not in resources:
                errors.append(f"resources must specify {var}")
            elif not isinstance(resources[var], str):
                errors.append(f"resources {var} must be an sequence (e.g. ACGT...)")
            elif len(set(resources[var]).difference("ACGT")) > 0:
                errors.append(f"resources {var} must contain only the letters ACGT")
            elif len(resources[var]) == 0:
                errors.append(f"resources {var} must contain at least one base")

        # Validate the sample barcodes
        for barcode in ['i5', 'i7']:
            missing = False
            for sample, vals in global_config['samples'].items():
                if barcode not in vals['barcodes']:
                    errors.append(f"samples.barcodes must specify {barcode} barcodes for sample {sample}")
                    missing = True
                elif not isinstance(vals['barcodes'][barcode], str):
                    errors.append(f"sample {sample} barcode {barcode} must be an sequence (e.g. ACGT...)")
                elif len(set(vals['barcodes'][barcode]).difference("ACGT")) > 0:
                    errors.append(f"sample {sample} barcode {barcode} must contain only the letters ACGT")
                elif len(vals['barcodes'][barcode]) == 0:
                    errors.append(f"sample {sample} barcode {barcode} must contain at least one base")
            if not missing:
                barcode_lengths = [len(sample_vals['barcodes'][barcode]) for sample_vals in global_config['samples'].values()]
                if len(set(barcode_lengths))  != 1:
                    errors.append(f"All samples {barcode} barcodes must be the same length.")

        return errors

    def read_flourescence(self, cluster: Cluster, range_start: int, range_end: int, direction: str, rng: np.random.Generator) -> tuple[np.ndarray, int, str, str]:
        """
        Generates flourescence readings from sequence-by-synthesis over the range given
        (specified from 5' to 3' if direction = '+' else from 3' to 5').
        Use call_bases to then get bases and quality scores from the read

        Parameters
        ----------

        cluster:
            cluster to read by sequence-by-synthesis
        range_start:
            The starting position (from the 5' end), 0-based
        range_end:
            The ending position (exlusive, from the 5' end), 0-based
        direction:
            '+' if reading from the 5'-3' direction, '-' if reading from 3'-5' direction
        rng:
            The random number generator

        Returns
        -------
            (flourescence, read_start, read_cigar, read_strand)
        """
        read_len = range_end - range_start

        assert cluster.base_counts is not None
        base_counts = cluster.base_counts
        if direction == '-':
            # Take reverse...
            base_counts = base_counts[:, ::-1]
            # ... and complement
            base_counts = base_counts[[3,2,1,0],:] # A<->T, C<->G

        # Grab the bases we care about
        # we go past the end of the read due to the possibility of skips (prephasing)
        base_counts = base_counts[:, range_start:range_end + self.MAX_SKIPS]


        # Pad with zeros for the possibility of skips/drops
        # that get past the start/end of the molecule
        end_pad_len = read_len + self.MAX_SKIPS - len(base_counts)
        prepad_len = self.MAX_DROPS
        padded_base_counts = np.concatenate(
                (   np.zeros((4,prepad_len)),
                    base_counts,
                    np.zeros((4,end_pad_len))
                ),
                axis=1)

        #TODO: this should all be done at the level of a lane or tile or something
        #      various parameters need to be estimated across all clusters

        # Add skips (aka prephasing) where bases are added too quickly in some strands
        # This causes bases to be read too early
        # We assume that each molecule skips at most MAX_SKIPS times
        # And we approximate skips and drops without simulating every molecule
        # individually by just using average molecule counts

        # Each base is equally likely to give a skip
        frac_skipped = get_frac_skipped_cython(self.skip_rate, self.MAX_SKIPS, cluster.molecule_count, read_len, rng)
        # Drops (aka (post)phasing), where bases are not added when they should have been,
        # are done the same as skips
        frac_dropped = get_frac_skipped_cython(self.drop_rate, self.MAX_DROPS, cluster.molecule_count, read_len, rng)

        # 'smear' base counts according to the number of skips
        smeared_base_counts = np.sum( [padded_base_counts[:, prepad_len + skip : prepad_len + read_len + skip] * frac_skipped[:,skip]
                                        for skip in range(1, self.MAX_SKIPS+1)],
                                      axis=0)
        # ...and by the number of drops
        smeared_base_counts = np.sum( [padded_base_counts[:, prepad_len - drop : prepad_len + read_len - drop] * frac_dropped[:,drop]
                                        for drop in range(1, prepad_len+1)],
                                      axis=0)
        # and the ones that neither drop nor skip
        frac_maintained = 1 - (1 - frac_skipped[:, 0]) - (1 - frac_dropped[:,0]) # Neither dropped nor skipped
        smeared_base_counts = padded_base_counts[:, prepad_len:prepad_len + read_len] * frac_maintained

        # Flourescence comes from these after including cross talk between the different colors
        flourescence = self.CROSS_TALK @ (smeared_base_counts + self.EPSILON * rng.normal(size=smeared_base_counts.shape))

        if direction == "+":
            # 1-based alignement
            align_start = range_start + 1
        else:
            align_start = len(cluster.molecule.sequence) - range_end

        read_start, read_cigar, read_strand = beers_utils.cigar.chain(
            align_start,
            f"{read_len}M",
            direction,
            cluster.molecule.source_start,
            cluster.molecule.source_cigar,
            cluster.molecule.source_strand,
        )
        return flourescence, read_start, read_cigar, read_strand

    def call_bases(self, flourescence: np.ndarray, epsilon_est: float, cross_talk_est_inv: np.ndarray) -> tuple[str, str]:
        '''
        From a flourescence reading, call sequences bases and quality score
        From an approximation of the Bustard algorithm

        Parameters
        ---------

        flourescence:
            2d array of shape (4, read_length) with flourescence values
            for each of the 4 frequencies for each base read
        epsilon_est:
            estimate for EPSILON, the noise size in the flourescence imaging
        cross_talk_est_inv:
            4x4 inverse of the estimate of the 4x4 cross talk matrix

        Returns
        -------
            called sequence and quality scores (as phred score string)
        '''

        ## We will approximate the Bustard algorithm for calling bases and quality scores
        # see https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC2765266/ for a description
        # of Bustard and some of the notation here.
        read_len = flourescence.shape[1]
        # gives probability of a template terminating at position j after t cycles
        phasing_matrix_inv = get_inv_phasing_matrix(read_len, self.skip_rate, self.drop_rate)
        base_counts_est = cross_talk_est_inv @ flourescence @ phasing_matrix_inv
        base_orders = np.argsort(base_counts_est, axis=0)
        called_base = base_orders[-1,:] # bases as nums 0,1,2,3 = ACGT
        second_best_base = base_orders[-2,:] # bases as nums 0,1,2,3 = ACGT
        M = (cross_talk_est_inv * epsilon_est) * (cross_talk_est_inv * epsilon_est).T
        highest_base_count = base_counts_est[called_base, np.arange(read_len)]
        second_base_count = base_counts_est[second_best_base, np.arange(read_len)]
        diff = highest_base_count - second_base_count
        contrasts = np.zeros(shape=(read_len, 4, 1))
        contrasts[np.arange(read_len), called_base] = 1
        contrasts[np.arange(read_len), second_best_base] = -1
        sigma = np.sqrt(contrasts.transpose(0,2,1) @ M @ contrasts)
        prob = scipy.stats.norm.sf(np.abs(diff) / sigma.flatten()) * 2 # two-tailed

        score = (-10 * np.log10(prob)).astype("uint8")
        score = np.minimum(PHRED_MAX_QUALITY, score)
        qual = (PHRED_MIN_ASCII + score).tobytes().decode()
        seq = BASE_ARRAY[called_base].tobytes().decode() # as string "ACGT..." 

        return seq, qual

@functools.lru_cache(maxsize=None)
def get_inv_phasing_matrix(read_len: int, skip_rate: float, drop_rate: float) -> np.ndarray:
    ''' Computes the inverse of Q, the phasing matrix

    The (j,t) entry of Q gives the probability of a template
    terminating at position j after t cycles

    Cached for speed-ups. Since reads are all the same length, they all get the same
    matrix and caching is very useful.

    Parameters
    ----------
    read_len:
        lenght of reads
    skip_rate:
        Rate of skipping events (forward phasing)
    drop_rate:
        Rate of skipping events (reverse phasing)

    Returns
    -------
        Inverse of the phasing matrix
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

def get_frac_skipped_py(rate: float, max_skips: int, molecule_count: int, read_len: int, rng: np.random.Generator) -> np.ndarray:
    ''' 
    Compute the array of length (read_len, max_skips+1) that indicate how many
    of the molecule_count molecules have incurred exactly i skips by the jth
    base in (i,j) element

    Parameters
    ----------
    rate:
        Positive number indicating rate of skipping
    max_skips:
        Maximum number of skipps we consider
    molecule_count:
        Number of molecules in the cluster
    read_len:
        lenght of reads
    rng:
        Random number generator instance
    '''
    if rate == 0:
        return np.zeros((read_len, max_skips+1))

    #num_skipped = np.zeros( (read_len, max_skips), dtype=int)
    #num_skipped[0,0] = molecule_count # Start with everything in no skips
    #for i in range(0, read_len):
    #    # Compute how many skips occur among the molecules
    #    # divided up by how many skips have already occured
    #    new_skips = rng.binomial(n=num_skipped[i, :-1], p=rate)
    #    num_skipped[i, 1:] += new_skips
    #    num_skipped[i, :-1] -= new_skips
    #    if i + 1 < read_len:
    #        num_skipped[i+1,:] = num_skipped[i,:]
    #frac_skipped = num_skipped / molecule_count


    # skip_locs[i,j] is the base number at which the ith molecule makes its jth skip
    # if past the end of read_len, then 'never' skips
    skip_locs = np.cumsum(rng.geometric(rate, size=(molecule_count, max_skips)), axis=1)
    # num_skipped[k,j] is the number of molecules that have had exactly j skips by the kth base
    num_skipped = np.count_nonzero(skip_locs[:,None,:] <= np.arange(1, read_len+1)[None,:,None], axis=0)
    frac_skipped = num_skipped / molecule_count
    # Add on the no-skipped ones
    frac_skipped = np.hstack([
        (np.ones(frac_skipped.shape[0]) * (1 - frac_skipped.sum(axis=1)))[:,None],
        frac_skipped
    ])
    #print(frac_skipped.shape, frac_skipped.dtype)
    #return frac_skipped


    #skips = rng.random(size = (molecule_count, read_len)) < rate
    #num_skips_so_far = np.minimum(np.cumsum(skips, axis=1), max_skips)
    ## Count number of molecules that have had exactly k skips at the nth base, for k= 0,1,...MAX_SKIPS
    #num_mols_skipped = np.count_nonzero((num_skips_so_far[:,:,None] == np.arange(0, max_skips+1)[None,None,:]), axis=0)
    #frac_skipped = num_mols_skipped / molecule_count
    #lkajf
    return frac_skipped
