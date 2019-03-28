#!/usr/bin/env python
import argparse
from camparee.camparee_controller import CampareeController

controller = CampareeController()

parser = argparse.ArgumentParser(description='CAMPAREE - RNA molecule simulator')
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument('-c', '--config', required=True, help='Full path to configuration file.')
optional_named = parser.add_argument_group('optional named arguments - these override configuration file arguments.')
optional_named.add_argument('-r', '--run_id', type=int, help="Integer used to specify run id.")
optional_named.add_argument('-d', '--debug', action='store_true',
                            help='Indicates whether additional diagnostics are printed.')
optional_named.add_argument('-m', '--dispatcher_mode', choices=['serial', 'multicore', 'lsf'],
                            help='Indicates whether to dispatch jobs serially, using multicore, or using lsf')

parser.set_defaults(func=controller.run_camparee_pipeline)

args = parser.parse_args()
args.func(args)
