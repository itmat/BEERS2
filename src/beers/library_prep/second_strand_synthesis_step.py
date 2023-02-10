import numpy as np

from beers_utils.molecule import Molecule
from beers_utils.general_utils import GeneralUtils
import beers_utils.cigar

class SecondStrandSynthesisStep:
    '''
    Second Strand Synthesis Step

    A library preparation step that places random primers onto the first strand cDNA
    and then extends the primers to obtain double-stranded DNA.

    Primer location can be biased by read sequence (see the 'position probability matrix' parameter).
    Even when not biased, priming can be either perfect or imperfect. Perfect means that the first
    base will always have a primer bind to it and therefore the entire molecule will be perfectly
    duplicated as cDNA. Imperfect means that primers are placed randomly (at a rate of 'primers_per_kb')
    and the 5'-most one is extended to create the molecule, potentially losing some of the 5' end of the
    molecule.
    '''

    name = "Second Strand Synthsis Step"

    def __init__(self, parameters, global_config):
        self.parameters = parameters
        self.global_config = global_config
        self.position_probability_matrix = self.get_position_probability_matrix(self.parameters)
        self.primer_length = self.position_probability_matrix.shape[1]
        self.primes_per_kb = self.parameters['primes_per_kb']
        self.perfect_priming = self.parameters['perfect_priming']
        print("Second Strand cDNA Synthesis Step instantiated.")

    def execute(self, molecule_packet, rng, log):
        print("Second strand cDNA synthesis step starting...")
        cdna_sample = []
        for molecule in molecule_packet.molecules:
            if len(molecule.sequence) < self.primer_length:
                continue # Can't prime a sequence this short, so it'll just drop

            if self.perfect_priming:
                # Perfectly prime at the first base, regardless of sequence
                primed_site = 0
            else:
                # First choose how many places will be primed
                number_of_primed_sites = rng.binomial(n = len(molecule.sequence), p = self.primes_per_kb/1000)
                if number_of_primed_sites == 0:
                    number_of_primed_sites = 1 # TODO: we force everything to be primed at least once, but maybe we should allow non-priming?

                # Then weight all possible priming spots according to how they match the given sequence
                seq_bases = GeneralUtils.sequence_to_matrix(molecule.sequence)
                weights = np.array([np.lib.stride_tricks.sliding_window_view(seq_bases[i], self.primer_length) * self.position_probability_matrix[i, :]
                                                        for i in range(4)]).sum(axis=0).prod(axis=1)

                if weights.sum() == 0:
                    # Should only occur if a sequence lacks all the normal bases, which is an error condition. We log it here
                    # and choose a random priming site
                    print(f"Weighting failed unexpectedly on the following molecule {molecule.molecule_id} with sequence {molecule.sequence}")
                    p = None
                else:
                    p = weights/weights.sum()

                # Then choose the priming sites and take the 5'-most one
                priming_sites = rng.choice(len(weights), p=p, size=number_of_primed_sites)
                primed_site = min(priming_sites)

            cdna_start = primed_site + 1
            cdna_cigar = f"{len(molecule.sequence)-cdna_start+1}M"
            cdna_source_start, cdna_source_cigar, cdna_source_strand = beers_utils.cigar.chain(
                    cdna_start, cdna_cigar, "-",
                    molecule.source_start, molecule.source_cigar, molecule.source_strand
            )
            cdna_seq = GeneralUtils.create_complement_strand(molecule.sequence[cdna_start - 1:])
            cdna_molecule = Molecule(
                    molecule.molecule_id + '.cdna' ,
                    cdna_seq,
                    start = cdna_start,
                    cigar = cdna_cigar,
                    transcript_id = molecule.transcript_id,
                    source_start = cdna_source_start,
                    source_cigar = cdna_source_cigar,
                    source_strand = cdna_source_strand,
                    source_chrom = molecule.source_chrom,
            )
            cdna_sample.append(cdna_molecule)
            log.write(cdna_molecule)

        print("Second strand cDNA synthesis step complete.")
        molecule_packet.molecules = cdna_sample
        return molecule_packet

    @staticmethod
    def get_position_probability_matrix(parameters):
        position_probability_matrix = np.array([
            parameters["position_probability_matrix"][base]
                for base in GeneralUtils.BASE_ORDER
        ])
        return position_probability_matrix


    @staticmethod
    def validate(parameters, global_config):
        errors = []
        if 'position_probability_matrix' not in parameters:
            errors.append("Must contain 'position_probability_matrix'")
        else:
            ppm = parameters['position_probability_matrix']
            okay = True
            length = None
            for base in GeneralUtils.BASE_ORDER:
                if base not in ppm:
                    errors.append(f"Must specify {base} in 'position_probability_matrix'")
                    continue
                probs = ppm[base]
                if not isinstance(probs, list) or not all(isinstance(x, (float, int)) for x in probs):
                    errors.append(f"Position probability matrix entry for {base} must be list of numbers")
                    okay = False
                    continue

                if length is None:
                    length = len(probs)
                elif length != len(probs):
                    errors.append(f"All position probability matrix entries must have the same length: had {len(probs)} for {base} but expected {length}")
                    okay = False
            if okay:
                ppm =  SecondStrandSynthesisStep.get_position_probability_matrix(parameters)
                sums = ppm.sum(axis=0)
                if not np.isclose(sums, 1).all():
                    errors.append(f"Position probability matrix entries must sum to 1 across all the four bases at each position")
                if not (ppm > 0).all():
                    errors.append("Position probability matrix entries must all be positive")

        if 'primes_per_kb' not in parameters:
            errors.append("Must contain 'primes_per_kb' parameter")
        elif not isinstance(parameters['primes_per_kb'], (float, int)):
            errors.append("primes_per_kb must be a number")
        elif parameters['primes_per_kb'] <= 0:
            errors.append("primer_per_kb must be positive")

        if 'perfect_priming' not in parameters:
            errors.append("Must specify 'perfect_priming' (as true/false)")
        elif not isinstance(parameters['perfect_priming'], bool):
            errors.append(f"Expected 'perfect_priming' to be a bool but instead was {parameters['perfect_priming']}")

        return errors
