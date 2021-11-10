from beers_utils.molecule import Molecule
from beers.utilities.library_prep_utils import Utils
import beers_utils.cigar


class FirstStrandSynthesisStep:

    name = "First Strand Synthsis Step"

    def __init__(self, log_file, parameters, global_config):
        self.history_filename = log_file
        self.parameters = parameters
        self.global_config = global_config
        print("First Strand cDNA Synthesis Step instantiated.")

    def execute(self, molecule_packet):
        """



        """
        print("First strand cDNA synthesis step starting...")
        cdna_sample = []
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in molecule_packet.molecules:
                cdna_seq = Utils.create_complement_strand(molecule.sequence)

                for primer in molecule.bound_molecules:
                    start = len(molecule.sequence) - (primer.start + len(primer))
                    end = start + len(primer)
                    cdna_seq = cdna_seq[:start] + primer.sequence + cdna_seq[end:]
                    #cdna_seq[start:end]  = primer.sequence

                first_primer = molecule.bound_molecules[0]
                # TODO: right start? For now, everything is primeed, so it should be
                cdna_start = 1
                cdna_cigar = f"{len(cdna_seq)}M"
                cdna_source_start, cdna_source_cigar = beers_utils.cigar.chain(
                    cdna_start, cdna_cigar, molecule.source_start, molecule.source_cigar
                )
                cdna_molecule = Molecule(molecule.molecule_id + '.cdna' ,cdna_seq[cdna_start:],
                        start = cdna_start,
                        cigar = cdna_cigar,
                        transcript_id = molecule.transcript_id,
                        source_start = molecule.source_start,
                        source_cigar = molecule.source_cigar,
                        source_strand = '-' if molecule.source_strand == '+' else '-',
                        source_chrom = molecule.source_chrom,
                )
                cdna_sample.append(cdna_molecule)
                log_file.write(cdna_molecule.log_entry())

            print("First strand cDNA synthesis step complete.")
            molecule_packet.molecules = cdna_sample
            return molecule_packet



    def validate(self):
        return True
