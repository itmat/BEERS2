#!/usr/bin/env python

import argparse
from beers.expression.expression_pipeline import ExpressionPipeline

parser = argparse.ArgumentParser(description='Run Expression Pipeline Only')
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument('-c', '--config', required=True, help='Full path to configuration file')
args = parser.parse_args()
ExpressionPipeline.main(args.config)

'''
Example Call:
python run_expression_pipeline_only.py -c ../config/config.json
'''