#!/usr/bin/env python3
import argparse

from beers.library_prep.pipeline import Pipeline

parser = argparse.ArgumentParser(description="Run the full pipeline with specified config file.")
parser.add_argument("--config", nargs=1, default="config/config.json")

args = parser.parse_args()

lib_prep = Pipeline(config_filename="config/config.json")

lib_prep.validate()
lib_prep.execute()
