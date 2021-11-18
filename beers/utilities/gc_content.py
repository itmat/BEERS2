#!/usr/bin/env python
import numpy as np
from beers_utils.molecule_packet import MoleculePacket

def gc_content(molecule):
    seq = molecule.sequence
    gc = (seq.count('C') + seq.count('G') + seq.count('N')/2) / len(seq)
    return gc

def packet_gc_content(molecule_packet):
    '''
    Assess the GC content distribution of a molecule packet

    Returns gc content per every 2.5%, i.e. percent molecules in each GC content bin
    of 0-2.5%, 2.5-5%,  ... 97.5-100%
    '''

    gcs = [gc_content(mol) for mol in molecule_packet.molecules]

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
