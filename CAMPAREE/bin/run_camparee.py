#!/usr/bin/env python
import argparse
from camparee.camparee_controller import CampareeController
from beers_utils.constants import SUPPORTED_SCHEDULER_MODES
from camparee.camparee_utils import CAMPAREE_VERSION

controller = CampareeController()

parser = argparse.ArgumentParser(description=f"CAMPAREE - RNA molecule simulator (v{CAMPAREE_VERSION})")
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument('-c', '--config', required=True, help='Full path to configuration file.')
optional_named = parser.add_argument_group('optional named arguments - these override configuration file arguments.')
optional_named.add_argument('-r', '--run_id', help="Alphanumberic used to specify run id (letters, numbers and"
                                                   " underscores only).")
optional_named.add_argument('-d', '--debug', action='store_true',
                            help='Indicates whether additional diagnostics are printed.')
optional_named.add_argument('-m', '--scheduler_mode', choices=SUPPORTED_SCHEDULER_MODES,
                            help='Indicates whether to dispatch jobs serially, using multicore, or using lsf')
optional_named.add_argument('-s', '--seed', type=int,
                            help='Optional integer value used as a seed for random number generation.')

parser.set_defaults(func=controller.run_camparee_pipeline)

args = parser.parse_args()
args.func(args)
