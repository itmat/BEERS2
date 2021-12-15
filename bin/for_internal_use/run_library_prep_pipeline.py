#!/usr/bin/env python
import argparse
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline

parser = argparse.ArgumentParser(description='Library Prep Pipeline')
parser.add_argument('-s', '--seed', required=True, help="Seed", type=int)
parser.add_argument('-c', '--config', required=True, help='Configuration')
parser.add_argument('-C', '--config_file', required=True, help='Configuration File Path')
parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
parser.add_argument('-p', '--molecule_packet_filename', required = False, default=None, help="Serialized molecule packet filename (or none, if using distribution).")
parser.add_argument('-d', '--directory_structure', required=True, help="Structure of data and logs directories.")
parser.add_argument('--packet_id', required=True, help="ID of the packet to be used (either int or None if using id from the packet")
parser.add_argument('-D', '--distribution_directory', required = False, default=None, help="Directory of CAMPAREE distributions (or none, if using molecule packet).")
parser.add_argument('-N', '--num_molecules', required = False, default=10_000, help="Number of molecules to generate if using distributions.", type=int)
parser.add_argument('-S', '--sample_id', required = False, default=None, help="Id number of the sample (used only when generating from distributions)", type=int)

args = parser.parse_args()
LibraryPrepPipeline.main(
        args.seed,
        args.config,
        args.config_file,
        args.output_directory,
        args.directory_structure,
        args.molecule_packet_filename,
        args.packet_id,
        distribution_directory = args.distribution_directory,
        molecules_per_packet_from_distribution = args.num_molecules,
        sample_id = args.sample_id,
)
