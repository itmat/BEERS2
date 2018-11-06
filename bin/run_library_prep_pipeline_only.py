#!/usr/bin/env python
import argparse

from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline

parser = argparse.ArgumentParser(description="Run the Library Prep Pipeline Only with specified config file.")
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument("-c", "--config", required=True, help='Full path to configuration file')
args = parser.parse_args()
LibraryPrepPipeline.main(args.config)
