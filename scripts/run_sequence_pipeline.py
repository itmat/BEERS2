#!/usr/bin/env python
import argparse
from beers.sequence.sequence_pipeline import SequencePipeline

parser = argparse.ArgumentParser(description='Sequence Pipeline')
parser.add_argument('-s', '--seed', required=True, help="Seed")
parser.add_argument('-c', '--config', required=True, help='Configuration')
parser.add_argument('-C', '--config_file', required=True, help='Configuration File Path')
parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
parser.add_argument('-p', '--cluster_packet_path', required=True, help="Serialized cluster packet file path.")
parser.add_argument('-d', '--directory_structure', required=True, help="Structure of data and logs directories.")

args = parser.parse_args()
SequencePipeline.main(
        args.seed,
        args.config,
        args.config_file,
        args.output_directory,
        args.directory_structure,
        args.cluster_packet_path,
)
