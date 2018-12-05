#!/usr/bin/env python
import argparse
import sys
from beers.controller import Controller

controller = Controller()

parser = argparse.ArgumentParser(description='BEERS Simulator, Version 2.0')
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument('-c', '--config', required=True, help='Full path to configuration file.')
optional_named = parser.add_argument_group('optional named arguments - these override configuration file arguments.')
optional_named.add_argument('-r', '--run_id', type=int, help="Integer used to specify run id.")
optional_named.add_argument('-d', '--debug', action='store_true',
                            help='Indicates whether additional diagnostics are printed.')
optional_named.add_argument('-m', '--dispatcher_mode', choices=['serial', 'multicore', 'lsf'],
                            help='Indicates whether to dispatch jobs serially, using multicore, or using lsf')
subparsers = parser.add_subparsers(help='pipeline subcommand', dest="subcommand")
subparsers.required = True

parser_expression_pipeline = subparsers.add_parser('expression_pipeline', help='Run the expression pipeline only')
parser_expression_pipeline.set_defaults(func=controller.run_expression_pipeline)

parser_library_prep_pipeline = subparsers.add_parser('library_prep_pipeline', help='Run the library prep pipeline only')
parser_library_prep_pipeline.set_defaults(func=controller.run_library_prep_pipeline)

parser_sequence_pipeline = subparsers.add_parser('sequence_pipeline', help='Run the sequence pipeline only')
parser_sequence_pipeline.set_defaults(func=controller.run_sequence_pipeline)

parser_prep_and_sequence_pipeline = subparsers.add_parser('prep_and_sequence_pipeline',
                                                          help='Run both the library prep and sequence pipelines')
parser_prep_and_sequence_pipeline.set_defaults(func=controller.run_prep_and_sequence_pipeline)

args = parser.parse_args()
args.func(args)
