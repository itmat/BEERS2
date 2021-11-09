#!/usr/bin/env python
import argparse
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline

parser = argparse.ArgumentParser(description='Library Prep Pipeline')
parser.add_argument('-s', '--seed', required=True, help="Seed")
parser.add_argument('-c', '--config', required=True, help='Configuration')
parser.add_argument('-C', '--config_file', required=True, help='Configuration File Path')
parser.add_argument('-i', '--input_directory', required=True, help='Path to input directory.')
parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
parser.add_argument('-p', '--molecule_packet_filename', required=True, help="Serialized molecule packet filename.")
parser.add_argument('-d', '--directory_structure', required=True, help="Structure of data and logs directories.")
parser.add_argument('--packet_id', required=True, help="ID of the packet to be used (either int or None if using id from the packet")

args = parser.parse_args()
LibraryPrepPipeline.main(
        args.seed,
        args.config,
        args.config_file,
        args.input_directory,
        args.output_directory,
        args.directory_structure,
        args.molecule_packet_filename,
        args.packet_id,
)
