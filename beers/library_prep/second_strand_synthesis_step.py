import numpy as np

from beers_utils.molecule import Molecule
from beers_utils.general_utils import GeneralUtils
import beers_utils.cigar


class SecondStrandSynthesisStep:

    name = "Second Strand Synthesis Step"

    def __init__(self, log_file, parameters, global_config):
        self.history_filename = log_file
        self.parameters = parameters
        self.global_config = global_config
        self.position_probability_matrix = np.array([
            self.parameters["position_probability_matrix"][base]
                for base in GeneralUtils.BASE_ORDER
        ])
        self.primer_length = self.position_probability_matrix.shape[1]
        self.primes_per_kb = self.parameters['primes_per_kb']
        self.perfect_priming = self.parameters['perfect_priming']
        print("Second Strand cDNA Synthesis Step instantiated.")

    def execute(self, molecule_packet):
        print("Second strand cDNA synthesis step starting...")
        # TODO: this takes ONLY the second synthesis strand. Is that correct?
        #       Is this right for non-strand-specific analyses?
        cdna_sample = []
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                # TODO: this code is duplicated with first strand synthesis - should we refactor it out?
                if len(molecule.sequence) < self.primer_length:
                    continue # Can't prime a sequence this short, so it'll just drop

                if self.perfect_priming:
                    # Perfectly prime at the first base, regardless of sequence
                    primed_site = 0
                else:
                    # First choose how many places will be primed
                    number_of_primed_sites = np.random.binomial(n = len(molecule.sequence), p = self.primes_per_kb/1000)
                    if number_of_primed_sites == 0:
                        number_of_primed_sites = 1 # TODO: we force everything to be primed at least once, but maybe we should allow non-priming?

                    # Then weight all possible priming spots according to how they match the given sequence
                    seq_bases = GeneralUtils.sequence_to_matrix(molecule.sequence)
                    weights = np.array([np.lib.stride_tricks.sliding_window_view(seq_bases[i], self.primer_length) * self.position_probability_matrix[i, :]
                                                            for i in range(4)]).sum(axis=0).prod(axis=1)
                    # Then choose the priming sites and take the 5'-most one
                    priming_sites = np.random.choice(len(weights), p=weights/weights.sum(), size=number_of_primed_sites)
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
                log_file.write(cdna_molecule.log_entry())

            print("Second strand cDNA synthesis step complete.")
            molecule_packet.molecules = cdna_sample
            return molecule_packet



    def validate(self):
        return True
