#!/usr/bin/env python
import numpy as np
from beers_utils.molecule_packet import MoleculePacket
from beers_utils import cigar

def gc_content(molecule, aligned_only):
    seq = molecule.sequence
    if aligned_only:
        # Extract just the parts of seq that aligned to reference
        i = 0
        segments = []
        for op, num in cigar.split_cigar(molecule.source_cigar):
            if cigar.consumes[op]['query']:
                if cigar.consumes[op]['ref']:
                    segments.append(seq[i:i+num])
                i += num
        seq = ''.join(segments)

    if len(seq) == 0:
        return float("NaN")
    gc = (seq.count('C') + seq.count('G') + seq.count('N')/2) / len(seq)
    return gc

def packet_gc_content(molecule_packet, aligned_only=False):
    '''
    Assess the GC content distribution of a molecule packet

    Returns gc content per every 2.5%, i.e. percent molecules in each GC content bin
    of 0-2.5%, 2.5-5%,  ... 97.5-100%

    aligned_only: if True, then restrict to the bases of the molecule sequence that
                    aligned to the reference sequence, according to the source_cigar
                    So, for example, polyA tails and adapters won't be counted.
                    Default is to count the entire molecule sequence
    '''

    gcs = [gc_content(mol, aligned_only) for mol in molecule_packet.molecules]

    bins = np.linspace(0, 1, 41)
    bins[-1] = 1.001 # make it right-inclusive for the last bin
    distribution,_ = np.histogram(gcs, bins)
    distribution = distribution / len(gcs)

    return distribution, bins

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Compute GC content of a molecule packet")
    parser.add_argument("molecule_packet")

    args = parser.parse_args()

    molecule_packet_file = args.molecule_packet
    molecule_packet = MoleculePacket.deserialize(molecule_packet_file)

    distribution, bins = packet_gc_content(molecule_packet)
    print('\n'.join(f"{bin:0.2%}\t{x}" for bin, x in zip(bins, distribution)))
