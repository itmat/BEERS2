#!/usr/bin/env python
import argparse
from beers.sequence.sequence_pipeline import SequencePipeline
import sys

parser = argparse.ArgumentParser(description='Sequence Pipeline')
parser.add_argument('-c', '--config', required=True, help='Configuration')
parser.add_argument('-i', '--input_directory', required=True, help='Path to input directory.')
parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
parser.add_argument('-p', '--cluster_packet_filename', required=True, help="Serialized Cluster Packet Filename.")
args = parser.parse_args()
SequencePipeline.main(args.config, args.input_directory, args.output_directory, args.cluster_packet_filename)