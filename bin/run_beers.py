#!/usr/bin/env python
import argparse
from beers.expression.expression_pipeline import ExpressionPipeline
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline


def run_expression_pipeline(args):
    ExpressionPipeline.main(args.config)

def run_library_prep_pipeline(args):
    LibraryPrepPipeline.main(args.config)


parser = argparse.ArgumentParser(description='BEERS Simulator, Version 2.0')
parser.add_argument('-c', '--config', required=True, help='Full path to configuration file')
subparsers = parser.add_subparsers(help='pipeline subcommand')

parser_expression_pipeline = subparsers.add_parser('expression_pipeline', help='Run the expression pipeline only')
parser_expression_pipeline.set_defaults(func=run_expression_pipeline)

parser_library_prep_pipeline = subparsers.add_parser('library_prep_pipeline', help='Run the library prep pipeline only')
parser_library_prep_pipeline.set_defaults(func=run_library_prep_pipeline)

args = parser.parse_args()
args.func(args)
