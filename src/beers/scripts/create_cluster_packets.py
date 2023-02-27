'''
Create cluster packets from molecule packets
'''

import resource
import json
import pathlib

import numpy as np

from beers_utils.molecule_packet import MoleculePacket
from beers.flowcell import Flowcell

configuration = json.loads(snakemake.params.configuration)

cluster_packet_paths = snakemake.output
molecule_packet_paths = snakemake.input

seed_list = [snakemake.params.seed, snakemake.wildcards.sample, snakemake.wildcards.packet_num, 3]
rng = np.random.default_rng(seed_list)

flowcell = Flowcell(configuration, rng)
valid, msg = flowcell.validate()
if not valid:
    raise (ControllerValidationException(msg))

# Generate the cluster packets
for molecule_packet_path, cluster_packet_path in zip(molecule_packet_paths, cluster_packet_paths):
    print(f"Loading {molecule_packet_path} to flowcell")
    molecule_packet = MoleculePacket.deserialize(molecule_packet_path)
    cluster_packet = flowcell.load_flowcell(molecule_packet)
    cluster_packet.serialize(cluster_packet_path)
    print(f"Output cluster packet to {cluster_packet_path}")
    print(f"Intermediate loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")
