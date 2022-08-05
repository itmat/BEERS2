'''
Create cluster packets from molecule packets
'''

import resource
import json
import pathlib

from beers_utils.molecule_packet import MoleculePacket
from beers.flowcell import Flowcell

stage_name = "sequence_pipeline"
configuration = json.loads(snakemake.params.configuration)

output_directory = pathlib.Path(snakemake.params.outdir)
output_directory.mkdir(exist_ok=True)

molecule_packet_file_paths = snakemake.input.packet_files_from_molecule_files + snakemake.input.packet_files_from_distribution

flowcell = Flowcell(configuration, configuration['flowcell'])
valid, msg = flowcell.validate()
if not valid:
    raise (ControllerValidationException(msg))

# Generate the cluster packets
cluster_packet_file_paths = []
for molecule_packet_filename in molecule_packet_file_paths:
    print(f"Loading {molecule_packet_filename} to flowcell")
    molecule_packet = MoleculePacket.deserialize(molecule_packet_filename)
    sample = molecule_packet.sample.sample_id
    sample_dir = output_directory / f"sample{sample}" 
    sample_dir.mkdir(exist_ok=True)

    cluster_packet = flowcell.load_flowcell(molecule_packet)
    cluster_packet_filename = f"cluster_packet_start_pkt{cluster_packet.cluster_packet_id}.gzip"
    cluster_packet_file_path =  sample_dir / cluster_packet_filename
    cluster_packet.serialize(cluster_packet_file_path)
    cluster_packet_file_paths.append(cluster_packet_file_path)
    print(f"Output cluster packet to {cluster_packet_file_path}")
    print(f"Intermediate loaded - process RAM at {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")
