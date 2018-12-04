#!/usr/bin/env python
import argparse
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
import sys

parser = argparse.ArgumentParser(description='Library Prep Pipeline')
parser.add_argument('-c', '--config', required=True, help='Configuration')
parser.add_argument('-i', '--input_directory', required=True, help='Path to input directory.')
parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
parser.add_argument('-p', '--molecule_packet_filename', required=True, help="Serialized Molecule Packet Filename.")
args = parser.parse_args()
LibraryPrepPipeline.main(args.config, args.input_directory, args.output_directory, args.molecule_packet_filename)