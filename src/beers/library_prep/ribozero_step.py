import sys
import re
import functools
from timeit import default_timer as timer
import numpy as np
import scipy.fft
from beers_utils.molecule import Molecule
from beers_utils.general_utils import GeneralUtils

class RiboZeroStep:
    """
    This step simulates the RiboZero procedure, which uses a library of oligos to deplete the ribosomal rRNA
    content of the sample.

    The library of oligos is sourced from:
        Adiconis, X., Borges-Rivera, D., Satija, R. et al.
        Comparative analysis of RNA sequencing methods for degraded or low-input samples.
        Nat Methods 10, 623-629 (2013). https://doi.org/10.1038/nmeth.2483

    Configuration Example::

        # Maximum chance to degrade a molecule that matches one of the oligos
        # (Molecules with multiple matching sites will still exceed this limit)
        max_degrade_chance: 0.99

        # Exponent to raise the score to to determine the chance of degradation
        # Exponent of 1 means that a 50% match has a 50% probability of degraded
        # exponent of 2 means a 50% match has a 25% probability to degrade, etc.
        # Note that 50% matches happen regularly and a 25% match happens ubiquitously
        # and this probability is per-base, so the exponent should not be low
        # E.g., a random length 2000 bp transcript will degrade at about 2.5 sites
        # at degrade_exponent=8. If degrade_entire_molecule=True, then the exponent
        # must be even higher (such as 15) to prevent losing most molecules since
        # any degradation event will remove the entire molecule.
        degrade_exponent: 8


        # If true, degrade/remove the entire molecule if any part matches
        # otherwise, fragments the molecule around the match and removes
        # the matching segment. Note that if this is True, then you must set the
        # degrade_exponent variable to be higher to prevent the loss of most molecules.
        degrade_entire_molecule: False
    """

    name = "RiboZero Step"

    def __init__(self, parameters, global_config):
        self.max_degrade_chance = parameters["max_degrade_chance"]
        self.degrade_exponent = parameters["degrade_exponent"]
        self.degrade_entire_molecule = parameters["degrade_entire_molecule"]
        self.global_config = global_config
        print("RiboZero selection step instantiated")

    def execute(self, molecule_packet, rng, log):
        print("RiboZero selection step starting")
        retained_molecules = []
        for molecule in molecule_packet.molecules:
            # Check if any of the oligo sequences are present in the molecule
            scores = ungapped_alignment_scores(ENCODED_OLIGO_LIBRARY, LIBRRAY_INDEX, molecule.sequence)
            degrade_probability = self.max_degrade_chance * (scores / OLIGO_LENGTH) ** self.degrade_exponent
            to_degrade = rng.random(size=degrade_probability.shape) <= degrade_probability
            degrade_sites, = np.where(to_degrade)

            times_degraded = 0
            remaining_molecule = molecule
            last_end = 1
            for site in degrade_sites:
                if site < last_end:
                    continue

                times_degraded += 1

                if self.degrade_entire_molecule:
                    # Don't retain any of the molecule
                    continue

                # Retain the piece between the end of the last fragment and this oligo match
                new_mol = molecule.make_fragment(last_end, site)
                if len(new_mol.sequence) > 0:
                    retained_molecules.append(new_mol)

                last_end = site + OLIGO_LENGTH
            if times_degraded > 0:
                note = f'degraded {times_degraded} times'

                if not self.degrade_entire_molecule:
                    # Take the piece after the last oligo match
                    new_mol = molecule.make_fragment(last_end, len(molecule.sequence))
                    if len(new_mol.sequence) > 0:
                        retained_molecules.append(new_mol)
            else:
                retained_molecules.append(molecule)
                note = 'retained'
            log.write(molecule, note)
        print("RiboZero selection step complete")
        molecule_packet.molecules = retained_molecules
        return molecule_packet

    @staticmethod
    def validate(parameters, global_config):
        errors = []
        if "max_degrade_chance" in parameters:
            max_degrade_chance = parameters["max_degrade_chance"]
            if max_degrade_chance < 0 or max_degrade_chance > 1:
                errors.append("The match degradation chance (max_degrade_chance) must be between 0 and 1")
        else:
            errors.append("Must specify max_degrade_chance (between 0 and 1)")

        if "degrade_exponent" in parameters:
            degrade_exponent = parameters["degrade_exponent"]
            if not (degrade_exponent > 0):
                errors.append("The degradation exponent (degrade_exponent) must be greater than 0")
        else:
            errors.append("Must specify degrade_exponent (a positive number)")

        if "degrade_entire_molecule" in parameters:
            degrade_entire_molecule = parameters["degrade_entire_molecule"]
            if degrade_entire_molecule not in [True, False]:
                errors.append("degrade_entire_molecule must be a boolean value")
        else:
            errors.append("Must specify degrade_entire_molecule (boolean)")

        return errors

ACGT = np.array([ord(x) for x in "ACGT"])
def encode(sequence):
    ''' Maps a sequence (as a string) to a numpy array of shape (4,n)
    which one-hot encodes the A/C/G/T value '''
    encoded = np.frombuffer(sequence.encode("ascii"), dtype='uint8')
    return (ACGT[:,None] == encoded[None,:]).astype(int)
def padded_encode(sequence, length):
    ''' Same as encode(sequence) but if sequence is less than given length
    it is padded with 1/1/1/1 values for ACGT to indicate it matches anything'''
    encoded = encode(sequence)
    return np.concatenate((
            encoded,
            np.full(shape=(4,max(0, length - len(sequence))), fill_value = 1),
        ),
        axis = 1,
    )

def ungapped_alignment_scores(queries, query_index, ref):
    ''' Returns the scores at each position of the best ungaped alignment
    of any of the query seqs against the reference sequence '''

    ref = encode(ref)
    scores = indexed_convolve(ref, queries, query_index)
    return scores.max(axis=(0,1))

def indexed_convolve(ref, library, index):
    '''
    Takes a library of queries and convolves them with ref

    The indexed library means an array of the fft's of
    all the library query values.
    The output is the (n, len(ref)-query_size) values of convolving
    each of the queries with the reference
    '''
    shape = [s1 + s2 - 1 for s1, s2 in zip(ref.shape, library.shape[1:])]
    shape_valid = [s1 - s2 + 1 for s1, s2 in zip(ref.shape, library.shape[1:])]

    # Get the index for the nearest 'fast' shape for the FFT
    # Using these shapes also reduces the number of precomputations
    # we have to do to form the library index since we use fewer different shapes
    fastshape = tuple(scipy.fft.next_fast_len(s) for s in shape)
    lib_index = index(fastshape)

    # Actually compute the convolution by multiplication in frequency space
    conv = scipy.fft.irfftn(
            scipy.fft.rfftn(ref, fastshape)[None, :, :] * lib_index,
            fastshape
    )

    # Cut out the part that we actually care about
    # End result will have the same shape as if you used method='valid'
    # in scipy.signal.convolve
    shape_slice = tuple(slice(sz) for sz in shape)
    ret = conv[(slice(None), *shape_slice)]
    return _centered(ret, (library.shape[0], *shape_valid))

def build_library_index(library):
    ''' For use with indexed_convolve.

    Computes fft of all the entires in a encoded sequence library
    '''
    @functools.cache
    def get_library_index(shape):
        return scipy.fft.rfftn(library[:, ::-1,::-1], shape, axes=(1,2))
    return get_library_index

def _centered(arr, newshape):
    # From scipy: https://github.com/scipy/scipy/blob/v1.10.1/scipy/signal/_signaltools.py#L1298-L1433
    # Return the center newshape portion of the array.
    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


OLIGO_LIBRARY = """
TAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACGACTTTTAC
TTCCTCTAGATAGTCAAGTTCGACCGTCTTCTCAGCGCTCCGCCAGGGCC
GTGGGCCGACCCCGGCGGGGCCGATCCGAGGGCCTCACTAAACCATCCAA
TCGGTAGTAGCGACGGGCGGTGTGTACAAAGGGCAGGGACTTAATCAACG
CAAGCTTATGACCCGCACTTACTCGGGAATTCCCTCGTTCATGGGGAATA
ATTGCAATCCCCGATCCCCATCACGAATGGGGTTCAACGGGTTACCCGCG
CCTGCCGGCGTAGGGTAGGCACACGCTGAGCCAGTCAGTGTAGCGCGCGT
GCAGCCCCGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCG
GGTGGCTGAACGCCACTTGTCCCTCTAAGAAGTTGGGGGACGCCGACCGC
TCGGGGGTCGCGTAACTAGTTAGCATGCCAGAGTCTCGTTCGTTATCGGA
ATTAACCAGACAAATCGCTCCACCAACTAAGAACGGCCATGCACCACCAC
CCACGGAATCGAGAAAGAGCTATCAATCTGTCAATCCTGTCCGTGTCCGG
GCCGGGTGAGGTTTCCCGTGTTGAGTCAAATTAAGCCGCAGGCTCCACTC
CTGGTGGTGCCCTTCCGTCAATTCCTTTAAGTTTCAGCTTTGCAACCATA
CTCCCCCCGGAACCCAAAGACTTTGGTTTCCCGGAAGCTGCCCGGCGGGT
CATGGGAATAACGCCGCCGCATCGCCGGTCGGCATCGTTTATGGTCGGAA
CTACGACGGTATCTGATCGTCTTCGAACCTCCGACTTTCGTTCTTGATTA
ATGAAAACATTCTTGGCAAATGCTTTCGCTCTGGTCCGTCTTGCGCCGGT
CCAAGAATTTCACCTCTAGCGGCGCAATACGAATGCCCCCGGCCGTCCCT
CTTAATCATGGCCTCAGTTCCGAAAACCAACAAAATAGAACCGCGGTCCT
ATTCCATTATTCCTAGCTGCGGTATCCAGGCGGCTCGGGCCTGCTTTGAA
CACTCTAATTTTTTCAAAGTAAACGCTTCGGGCCCCGCGGGACACTCAGC
TAAGAGCATCGAGGGGGCGCCGAGAGGCAAGGGGCGGGGACGGGCGGTGG
CTCGCCTCGCGGCGGACCGCCCGCCCGCTCCCAAGATCCAACTACGAGCT
TTTTAACTGCAGCAACTTTAATATACGCTATTGGAGCTGGAATTACCGCG
GCTGCTGGCACCAGACTTGCCCTCCAATGGATCCTCGTTAAAGGATTTAA
AGTGGACTCATTCCAATTACAGGGCCTCGAAAGAGTCCTGTATTGTTATT
TTTCGTCACTACCTCCCCGGGTCGGGAGTGGGTAATTTGCGCGCCTGCTG
CCTTCCTTGGATGTGGTAGCCGTTTCTCAGGCTCCCTCTCCGGAATCGAA
CCCTGATTCCCCGTCACCCGTGGTCACCATGGTAGGCACGGCGACTACCA
TCGAAAGTTGATAGGGCAGACGTTCGAATGGGTCGTCGCCGCCACGGG
GCGTGCGATCGGCCCGAGGTTATCTAGAGTCACCAAAGCCGCCGGCGCCC
GCCCCCCGGCCGGGGCCGGAGAGGGGCTGACCGGGTTGGTTTTGATCTGA
TAAATGCACGCATCCCCCCCGCGAAGGGGGTCAGCGCCCGTCGGCATGTA
TTAGCTCTAGAATTACCACAGTTATCCAAGTAGGAGAGGAGCGAGCGACC
AAAGGAACCATAACTGATTTAATGAGCCATTCGCAGTTTCACTGTACCGG
CCGTGCGTACTTAGACATGCATGGCTTAATCTTTGAGACAAGCATATGCT
TGGCTTAATCTTTGAGACAAGCATATGCTACTGGCAGGATCAACCAGGTA
GACAAACCCTTGTGTCGAGGGCTGACTTTCAATAGATCGCAGCGAGGGAG
CTGCTCTGCTACGTACGAAACCCCGACCCAGAAGCAGGTCGTCTACGAAT
GGTTTAGCGCCAGGTTCCCCACGAACGTGCGGTGCGTGACGGGCGAGGG
GCGGCCGCCTTTCCGGCCGCGCCCCGTTTCCCAGGACGAAGGGCACTCCG
CACCGGACCCCGGTCCCGGCGCGCGGCGGGGCACGCGCCCTCCCGCGGCG
GGGCGCGTGGAGGGGIGGGCGGCCCGCCGGCGGGGACAGGCGGGGGACCG
GCTATCCGAGGCCAACCGAGGCTCCGCGGCGCTGCCGTATCGTTCGCCTG
GGCGGGATTCTGACTTAGAGGCGTTCAGTCATAATCCCACAGATGGTAGC
TTCGCCCCATTGGCTCCTCAGCCAAGCACATACACCAAATGTCTGAACCT
GCGGTTCCTCTCGTACTGAGCAGGATTACCATGGCAACAACACATCATCA
GTAGGGTAAAACTAACCTGTCTCACGACGGTCTAAACCCAGCTCACGTTC
CCTATTAGTGGGTGAACAATCCAACGCTTGGCGAATTCTGCTTCACAATG
ATAGGAAGAGCCGACATCGAAGGATCAAAAAGCGACGTCGCTATGAACGC
TTGGCCGCCACAAGCCAGTTATCCCTGTGGTAACTTTTCTGACACCTCCT
GCTTAAAACCCAAAAGGTCAGAAGGATCGTGAGGCCCCGCTTTCACGGTC
TGTATTCGTACTGAAAATCAAGATCAAGCGAGCTTTTGCCCTTCTGCTCC
ACGGGAGGTTTCTGTCCTCCCTGAGCTCGCCTTAGGACACCTGCGTTACC
GTTTGACAGGTGTACCGCCCCAGTCAAACTCCCCACCTGGCACTGTCCCC
GGAGCGGGTCGCGCCCGGCCGGGCGGGCGCTTGGCGCCAGAAGCGAGAGC
CCCTCGGGCTCGCCCCCCCGCCTCACCGGGTCAGTGAAAAAACGATCAGA
GTAGTGGTATTTCACCGGCGGCCCGCAGGGCCGCGGACCCCGCCCCGGGC
CCCTCGCGGGGACACCGGGIGGGCGCCGGGGGCCTCCCACTTATTCTACA
CCTCTCATGTCTCTTCACCGTGCCAGACTAGAGTCAAGCTCAACAGGGTC
TTCTTTCCCCGCTGATTCCGCCAAGCCCGTTCCCTTGGCTGTGGTTTCGC
TGGATAGTAGGTAGGGACAGTGGGAATCTCGTTCATCCATTCATGCGCGT
CACTAATTAGATGACGAGGCATTTGGCTACCTTAAGAGAGTCATAGTTAC
TCCCGCCGTTTACCCGCGCTTCATTGAATTTCTTCACTTTGACATTCAGA
GCACTGGGCAGAAATCACATCGCGTCAACACCCGCCGCGGGCCTTCGCGA
TGCTTTGTTTTAATTAAACAGTCGGATTCCCCTGGTCCGCACCAGTTCTA
AGTCGGCTGCTAGGCGCCGGCCGAGGCGAGGCGCGCGCGGAACCGCGGCC
CCGGGGGCGGACCCGGCGGGIGGGACCGGCCCGCGGCCCCTCCGCCGCCT
GCCGCCGCCGCCGCCGCGCGCCGAGGAGGAGGGGGGAACGGGGGGCGGAC
GGGCCGGGIGGGTAGGGCGGGGGGACGAACCGCCCCGCCCCGCCGCCCG
CCGACCGCCGCCGCCCGACCGCTCCCGCCCCCAGCGGACGCGCGCGCGAC
CGAGACGTGGGGTGGGGGTGGGGGGCGCGCCGCGCCGCCGCCGGGCTCCC
CGGGGGCGGCCGCGACGCCCGCCGCAGCTGGGGCGATCCACGGGAAGGGC
CCGGCTCGCGTCCAGAGTCCGCGCCGCCGCCGGCCCCCCGGGTCCCCGGG
GCCCCCCTCGCGGGGACCTGCCCCCGCCGGCCGCCCCGGCGGCCGCCGCG
CGGCCCCTGCCGCCCCGACCCTTCTCCCCCCGCCGCGCCCCCACGCGGCG
CTCCCCCGGGGAGGGGGGAGGACGGGGAGCGGGGGAGAGAGAGAGAGAGA
GGGCGCGGGGTGGGGAGGGAGCGAGCGGCGCGCGCGGGTGGGGCGGGGGA
GGGCCGCGAGGGGGGTGCCCCGGGCGTGGGGIGGGCGCGCGCCTCGTCCA
GCCGCGGCGCGCGCCCAGCCCCGCTTCGCGCCCCAGCCCGACCGACCCAG
CCCTTAGAGCCAATCCTTATCCCGAAGTTACGGATCCGGCTTGCCGACTT
CCCTTACCTACATTGTTCCAACATGCCAGAGGCTGTTCACCTTGGAGACC
TGCTGCGGATATGGGTACGGCCCGGCGCGAGATTTACACCCTCTCCCCCG
GATTTTCAAGGGCCAGCGAGAGCTCACCGGACGCCGCCGGAACCGCGACG
CTTTCCAAGGCACGGGCCCCTCTCTCGGGGCGAACCCATTCCAGGGCGCC
CTGCCCTTCACAAAGAAAAGAGAACTCTCCCCGGGGCTCCCGCCGGCTTC
TCCGGGATCGGTCGCGTTACCGCACTGGACGCCTCGCGGCGCCCATCTCC
GCCACTCCGGATTCGGGGATCTGAACCCGACTCCCTTTCGATCGGCCGAG
GGCAACGGAGGCCATCGCCCGTCCCTTCGGAACGGCGCTCGCCCATCTCT
CAGGACCGACTGACCCATGTTCAACTGCTGTTCACATGGAACCCTTCTCC
ACTTCGGCCTTCAAAGTTCTCGTTTGAATATTTGCTACTACCACCAAGAT
CTGCACCTGCGGCGGCTCCACCCGGGCCCGCGCCCTAGGCTTCAAGGCTC
ACCGCAGCGGCCCTCCTACTCGTCGCGGCGTAGCGTCCGCGGGGCTCCGG
GGGCGGGGAGCGGGGCGTGGGCGGGAGGAGGGGAGGAGGCGTGGG
GGGCGGGGGAAGGACCCCACACCCCCGCCGCCGCCGCCGCCGCCGCCCTC
CGACGCACACCACACGCGCGCGCGCGCGCGCCGCCCCCGCCGCTCCCGTC
CACTCTCGACTGCCGGCGACGGCCGGGTATGGGCCCGACGCTCCAGCGCC
ATCCATTTTCAGGGCTAGTTGATTCGGCAGGTGAGTTGTTACACACTCCT
TAGCGGATTCCGACTTCCATGGCCACCGTCCTGCTGTCTATATCAACCAA
CACCTTTTCTGGGGTCTGATGAGCGTCGGCATCGGGCGCCTTAACCCGGC
GTTCGGTTCATCCCGCAGCGCCAGTTCTGCTTACCAAAAGTGGCCCACTA
GGCACTCGCATTCCACGCCCGGCTCCACGCCAGCGAGCCGGGCTTCTTAC
CCATTTAAAGTTTGAGAATAGGTTGAGATCGTTTCGGCCCCAAGACCTCT
AATCATTCGCTTTACCGGATAAAACTGCGTGGCGGGGGTGCGTCGGGTCT
GCGAGAGCGCCAGCTATCCTGAGGGAAACTTCGGAGGGAACCAGCTACTA
GATGGTTCGATTAGTCTTTCGCCCCTATACCCAGGTCGGACGACCGATTT
GCACGTCAGGACCGCTACGGACCTCCACCAGAGTTTCCTCTGGCTTCGCC
CTGCCCAGGCATAGTTCACCATCTTTCGGGTCCTAACACGTGCGCTCGTG
CTCCACCTCCCCGGCGCGGCGGGCGAGACGGGCCGGTGGTGCGCCCTCGG
CGGACTGGAGAGGCCTCGGGATCCCACCTCGGCCGGCGAGCGCGCCGGCC
TTCACCTTCATTGCGCCACGGCGGCTTTCGTGCGAGCCCCCGACTCGCGC
ACGTGTTAGACTCCTTGGTCCGTGTTTCAAGACGGGTCGGGTGGGTAGCC
GACGTCGCCGCCGACCCCGTGCGCTCGCTCCGCCGTCCCCCTCTTCGGG
GACGCGCGCGTGGCCCCGAGAGAACCTCCCCCGGGCCCGACGGCGCGACC
CGCCCGGGGCGCACTGGGGACAGTCCGCCCCGCCCCCCGACCCGCGCGCG
GCACCCCCCCCGTCGCCGGGGCGGGGGCGCGGGGAGGAGGGGTGGGAGAG
CGGTCGCGCCGTGGGAGGGGTGGCCCGGCCCCCCCACGAGGAGACGCCGG
CGCGCCCCCGCGGGGGAGACCCCCCTCGCGGGGGATTCCCCGCGGGGGTG
GGCGCCGGGAGGGGGGAGAGCGCGGCGACGGGTCTCGCTCCCTCGGCCCC
GGGATTCGGCGAGTGCTGCTGCCGGGGGGGCTGTAACACTCGGGGIGGGT
TTCGGTCCCGCCGCCCCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCGC
CCCGACCCGCGCGCCCTCCCGAGGGAGGACGCGGGGCCGGGGGGCGGAGA
CGGGGGAGGAGGAGGACGGACGGACGGACGGGGCCCCCCGAGCCACCTTC
CCCGCCGGGCCTTCCCAGCCGTCCCGGAGCCGGTCGCGGCGCACCGCCGC
GGTGGAAATGCGCCCGGCGGCGGCCGGTCGCCGGTCGGGGGACGGTCCCC
CGCCGACCCCACCCCCGGCCCCGCCCGCCCACCCCCGCACCCGCCGGAGC
CCGCCCCCTCCGGGGAGGAGGAGGAGGGGCGGCGGGGGAAGGGAGGGCGG
GTGGAGGGGTCGGGAGGAACGGGGGGCGGGAAAGATCCGCCGGGCCGCCG
ACACGGCCGGACCCGCCGCCGGGTTGAATCCTCCGGGCGGACTGCGCGGA
CCCCACCCGTTTACCTCTTAACGGTTTCACGCCCTCTTGAACTCTCTCTT
CAAAGTTCTTTTCAACTTTCCCTTACGGTACTTGTTGACTATCGGTCTCG
TGCCGGTATTTAGCCTTAGATGGAGTTTACCACCCGCTTTGGGCTGCATT
CCCAAGCAACCCGACTCCGGGAAGACCCGGGCGCGCGCCGGCCGCTACCG
GCCTCACACCGTCCACGGGCTGGGCCTCGATCAGAAGGACTTGGGCCCCC
CACGAGCGGCGCCGGGGAGCGGGTCTTCCGTACGCCACATGTCCCGCGCC
CCGCGGGGCGGGGATTCGGCGCTGGGCTCTTCCCTGTTCACTCGCCGTTA
CTGAGGGAATCCTGGTTAGTTTCTTTTCCTCCGCTGACTAATATGCTTAA
GACTAATATGCTTAAATTCAGCGGGTCGCCACGTCTGATCTGAGGTCGCG
AAGCGACGCTCAGACAGGCGTAGCCCCGGGAGGAACCCGGGGCCGCAAGT
GCGTTCGAAGTGTCGATGATCAATGTGTCCTGCAATTCACATTAATTCTC
GCAGCTAGCTGCGTTCTTCATCGACGCACGAGCCGAGTGATCCACCGCTA
AAACCCTGTTCTTGGGTGGGTGTGGGTATAATACTAAGTTGAGATGATAT
CATTTACGGGGGAAGGCGCTTTGTGAAGTAGGCCTTATTTCTCTTGTCCT
TTCGTACAGGGAGGAATTTGAANGTAGATAGAAACCGACCTGGATTACTC
CGGTCTGAACTCAGATCACGTAGGACTTTAATCGTTGAACAAACGAACCT
TTAATAGCGGCTGCACCATCGGGATGTCCTGATCCAACATCGAGGTCGTA
AACCCTATTGTTGATATGGACTCTAGAATAGGATTGCGCTGTTATCCCTA
GGGTAACTTGTTCCGTTGGTCAAGTTATTGGATCAATTGAGTATAGTAGT
TCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTGGGTTCT
GCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAG
GACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAATTAAAGCTCCAT
AGGGTCTTCTCGTCTTGCTGTGTTATGCCCGCCTCTTCACGGGCAGGTCA
ATTTCACTGGTTAAAAGTAAGAGACAGCTGAACCCTCGTGGAGCCATTCA
TACAGGTCCCTATTTAAGGAACAAGTGATTATGCTACCTTTGCACGGTTA
GGGTACCGCGGCCGTTAAACATGTGTCACTGGGCAGGCGGTGCCTCTAAT
ACTGGTGATGCTAGAGGTGATGTTTTTGGTAAACAGGCGGGGTAAGATTT
GCCGAGTTCCTTTTACTTTTTTTAACCTTTCCTTATGAGCATGCCTGTGT
TGGGTTGACAGTGAGGGTAATAATGACTTGTTGGTTGATTGTAGATATTG
GGCTGTTAATTGTCAGTTCAGTGTTTTAATCTGACGCAGGCTTATGCGGA
GGAGAATGTTTTCATGTTACTTATACTAACATTAGTTCTTCTATAGGGTG
ATAGATTGGTCCAATTGGGTGTGAGGAGTTCAGTTATATGTTTGGGATTT
TTTAGGTAGTGGGTGTTGAGCTTGAACGCTTTCTTAATTGGTGGCTGCTT
TTAGGCCTACTATGGGTGTTAAATTTTTTACTCTCTCTACAAGGTTTTTT
CCTAGTGTCCAAAGAGCTGTTCCTCTTTGGACTAACAGTTAAATTTACAA
GGGATTTAGAGGGTTCTGTGGGCAAATTTAAAGTTGAACTAAGATTCTA
TCTTGGACAACCAGCTATCACCAGGCTCGGTAGGTTTGTCGCCTCTACCT
ATAAATCTTCCCACTATTTTGCTACATAGACGGGTGTGCTCTTTTAGCTG
TTCTTAGGTAGCTCGTCTGGTTTCGGGGGTCTTAGCTTTGGCTCTCCTTG
CAAAGTTATTTCTAGTTAATTCATTATGCAGAAGGTATAGGGGTTAGTCC
TTGCTATATTATGCTTGGTTATAATTTTTCATCTTTCCCTTGCGGTACTA
TATCTATTGCGCCAGGTTTCAATTTCTATCGCCTATACTTTATTTGGGTA
AATGGTTTGGCTAAGGTTGTCTGGTAGTAAGGTGGAGTGGGTTTGGGGCT
GTTCGTCCAAGTGCACTTTCCAGTACACTTACCATGTTACGACTTGTCTC
CTCTATATAAATGCGTAGGGGTTTTAGTTAAATGTCCTTTGAAGTATACT
TGAGGAGGGTGACGGGCGGTGTGTACGCGCTTCAGGGCCCTGTTCAACTA
AGCACTCTACTCTTAGTTTACTGCTAAATCCACCTTCGACCCTTAAGTTT
CATAAGGGCTATCGTAGTTTTCTGGGGTAGAAAATGTAGCCCATTTCTTG
CCACCTCATGGGCTACACCTTGACCTAACGTCTTTACGTGGGTACTTGCG
CTTACTTTGTAGCCTTCATCAGGGTTTGCTGAAGATGGCGGTATATAGGC
TGAGCAAGAGGTGGTGAGGTTGATCGGGGTTTATCGATTACAGAACAGGC
TCCTCTAGAGGGATATGAAGCACCGCCAGGTCCTTTGAGTTTTAAGCTGT
GGCTCGTAGTGTTCTGGCGAGCAGTTTTGTTGATTTAACTGTTGAGGTTT
AGGGCTAAGCATAGTGGGGTATCTAATCCCAGTTTGGGTCTTAGCTATTG
TGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTG
GAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGAT
CTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTGAC
CGCGGTGGCTGGCACGAAATTGACCAACCCTGGGGTTAGTATAGCTTAGT
TAAACTTTCGTTTATTGCTAAAGGTTAATCACTGCTGTTTCCCGTGGG
TGTGGCTAGGCTAAGCGTTTTGAGCTGCATTGCTGCGTGCTTGATGCTTG
TTCCTTTTGATCGTGGTGATTTAGAGGGTGAACTCACTGGAACGGGGATG
CTTGCATGTGTAATCTTACTAAGAGCTAATAGAAAGGCTAGGACCAAACC
AAAGCCTACAGCACCCGGTATTCCCAGGCGGTCTCCCATCCAAGTACTAA
CCAGGCCCGACCCTGCTTAGCTTCCGAGATCAGACGAGATCGGGCGCGTT
TTCCGAGATCAGACGAGATCGGGCGCGTTCAGGGTGGTATGGCCGTAGAC
""".strip().split()
OLIGO_LENGTH = max(len(oligo) for oligo in OLIGO_LIBRARY)
REVERSE_OLIGO_LIBRARY = [GeneralUtils.create_complement_strand(oligo) for oligo in OLIGO_LIBRARY]
# We want to search for the reverse complement of the oligos in the molecules
ENCODED_OLIGO_LIBRARY = np.array([padded_encode(oligo, OLIGO_LENGTH) for oligo in REVERSE_OLIGO_LIBRARY])
LIBRRAY_INDEX = build_library_index(ENCODED_OLIGO_LIBRARY)