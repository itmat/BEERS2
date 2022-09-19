import argparse
import importlib
import pathlib
import time
import re
import os
import sys
import resource
import json
import math

import numpy as np

from beers_utils.sample import Sample
from beers_utils.constants import CONSTANTS
from beers_utils.molecule import Molecule
from beers_utils.molecule_packet import MoleculePacket
from beers_utils.general_utils import GeneralUtils

from camparee.molecule_maker import MoleculeMakerStep

from beers.abstract_beers_pipeline import AbstractBeersPipeline

class LibraryPrepPipeline():
    """
    The class runs all the steps in the library_prep pipeline as described and wired together in the configuration
    file.  The point of entry into this class is the static method main().  In HPC scenarios, this object's main is
    invoked on each node included in the compute, with a different molecule packet delivered to each node.
    """

    stage_name = "library_prep_pipeline"
    pipeline_log_subdirectory_name = "Pipeline"
    package = "beers.library_prep"

    def __init__(self):
        pass

    @staticmethod
    def validate(configuration, global_config):
        """
        Static method to run each step validate process to identify errant parameters.  If any errors are found,
        a validation exception is raised.
        """
        #if not all([molecule.validate() for molecule in self.molecule_packet.molecules]):
        #    raise BeersLibraryPrepValidationException("Validation error in molecule packet: see stderr for details.")

        # Gather the step classes
        validation_okay = True
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)

            # Valdiate steps's parameters and print out their errors (if any)
            errors = step_class.validate(parameters, global_config)
            if len(errors) > 0:
                validation_okay = False
                print(f"Validation errors in {step['step_name']}:", file=sys.stderr)
                for error in errors:
                    print('\t' + error, file=sys.stderr)

        if not validation_okay:
            raise BeersLibraryPrepValidationException("Validation error: see stderr for details.")

    def execute(self, configuration, global_config, output_directory, log_directory, input_molecule_packet, rng):
        """
        Opens the pipeline log for writing and serially runs the execute method of each step object found in the
        step list generated when this pipeline stage was initialized.  The final product (a modified molecule
        packet) is serialized into a data file.
        Note that the directory structure
        has already been created by the controller.  The steps described in the configuration dictionary are
        instantiated and those instantiated steps are added to a list for later invocation.
        :param configuration:  dictionary of the configuration data relevant to this pipeline stage.
        :param global_config: dictionary of full config data
        :param output_directory: file path to directory where results (molecule packet and quant file) should be output to
        :param log_directory: file path to write the log file to
        :param input_molecule_packet: the molecule packet to run through this pipeline stage
        :param rng: random number generator object
        """

        original_ids = set(str(m.molecule_id) for m in input_molecule_packet.molecules)
        self.print_summary(original_ids, input_molecule_packet.molecules)

        packet_num = input_molecule_packet.molecule_packet_id
        output_packet_path = pathlib.Path(output_directory) / f"library_prep_pipeline_result_molecule_pkt{packet_num}.txt"
        output_quant_path = pathlib.Path(output_directory) / f"library_prep_pipeline_result_molecule_pkt{packet_num}.quant_file"
        input_quant_file_path = pathlib.Path(log_directory) / f"library_prep_pipeline_molecule_pkt{packet_num}.quant_file"

        # Initialize steps prior to executing them below.
        pipeline_steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"molecule_pkt{input_molecule_packet.molecule_packet_id}.log"
            step_log_dir = pathlib.Path(log_directory) / step_name
            step_log_dir.mkdir(exist_ok=True)
            step_log_file_path = step_log_dir / step_log_filename
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)
            pipeline_steps.append(step_class(step_log_file_path, parameters, global_config))

        pipeline_start = time.time()

        molecule_packet = input_molecule_packet
        molecule_packet.write_quantification_file(input_quant_file_path)

        for step in pipeline_steps:
            step_name = step.__class__.name if hasattr(step.__class__, 'name') else step.__class__.__name__
            step_start = time.time()
            molecule_packet = step.execute(molecule_packet, rng)
            elapsed_time = time.time() - step_start

            self.print_summary(original_ids, molecule_packet.molecules, elapsed_time)

        pipeline_elapsed_time = time.time() - pipeline_start
        print(f"Finished {LibraryPrepPipeline.stage_name} in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a gzip file for inspection
        molecule_packet.serialize(output_packet_path)
        print(f"Output final sample to {output_packet_path}")
        molecule_packet.write_quantification_file(output_quant_path)
        print(f"Output final sample quantification to {output_quant_path}")

    def print_summary(self, original_ids, sample, elapsed_time=None):
        """Output a summary of the sample (number of molecules, time taken, etc.).

        :param sample: Molecule IDs from input molecule packet
        :param sample:
        :param elapsed_time:
        """
        if elapsed_time is not None:
            print(f"Step took {elapsed_time:.3} seconds")

        print(f"Sample has {len(sample)} molecules")
        parent_ids = set(str(m.molecule_id).split(".")[0] for m in sample)
        #TODO: this does not seem to be calculated correctly: always gives 0
        percent_original_represented = len(original_ids.intersection(parent_ids))/len(original_ids)
        print(f"Percent of the original ids that are still represented: {percent_original_represented:0.2%}")

        size_bin_cutoffs = [0, 100, 500, 1000, 1_000_000_000]
        size_counts, _ = np.histogram([len(mol) for mol in sample], bins=size_bin_cutoffs)
        print(f"Counts of molecules in size ranges:")
        for i in range(len(size_counts)-1):
            print(f" <={size_bin_cutoffs[i+1]}: {size_counts[i]}")
        print(f">{size_bin_cutoffs[-2]}: {size_counts[-1]}")

    @staticmethod
    def main(seed, configuration, global_configuration, output_directory, log_directory,
             molecule_packet_filename, packet_id, distribution_directory=None,
             molecules_per_packet_from_distribution=10000, sample_id=None,
             ):
        """
        This method would be called by a command line script in the bin directory.  It sets a random seed, loads a
        directory containing the relevant parts of the user's configuration file, unmarshalls a molecule packet from
        the provided molecule packet filename, initializes and validates the library prep pipeline stage and then
        executes it for the molecule packet.
        :param seed: value to use as the seed for the random number generator
        :param configuration: the json string containing the configration data specific to the library prep pipeline
        :param global_configuration: json string containing the configuration data for the whole run
        :param output_directory: path to directory where we will output molecule files and quant files
        :param log_directory: path to log directories
        :param directory_structure: instructions for creating the scaffolding needed to house the pipeline data and logs
        :param molecule_packet_filename: the file from which to unmarshall the molecule packet
        :param packet_id: id number to assign the packet
        :param distribution_directory: directory of CAMPAREE output distribution datas from which to generate molecules
            (default None, must supply molecule_packet_filename instead)
        :param molecules_per_packet_from_distribution: packet size to generate if using distributions
        :param sample_id: id number of the sample (if None, derive from the molecule packet provided)
        """

        configuration = json.loads(configuration)
        global_configuration = json.loads(global_configuration)
        molecule_maker_parameters = global_configuration['molecule_maker_parameters']
        packet_id = None if packet_id == 'None' else int(packet_id)

        # Seed the RNG by the given (global) seed, plus the sample ID, the packet ID, and a constant for the library prep stage of 1
        # This ensures that each packet is repeatable but that each gets their own seed state and that the seed state differs
        # from each of the stages it goes through.
        seed_list = [seed, sample_id, packet_id, 1]
        print(f"Initializing with seed {seed_list}")
        rng = np.random.default_rng(seed_list)

        if molecule_packet_filename:
            molecule_packet = MoleculePacket.from_CAMPAREE_molecule_file(molecule_packet_filename, packet_id)
        elif distribution_directory:
            distribution_dir = pathlib.Path(distribution_directory)
            sample = Sample(
                sample_id = sample_id,
                sample_name = sample_id,
                adapter_sequences = [],
                fastq_file_paths = [],
                pooled = False,
            )
            molecule_maker = MoleculeMakerStep(
                    log_directory_path = log_directory,
                    parameters = molecule_maker_parameters,
            )
            molecule_packets = list(molecule_maker.execute(
                    sample = sample,
                    sample_data_directory = distribution_dir,
                    output_type =  "generator",
                    output_molecule_count = molecules_per_packet_from_distribution, # we only generate one packet
                    seed = seed,
                    molecules_per_packet = molecules_per_packet_from_distribution,
                    rng = rng,
            ))
            molecule_packet = molecule_packets[0]
            molecule_packet.molecule_packet_id = packet_id
        else:
            raise BeersLibraryPrepValidationException("Neither molecule packet filename nor distribution directory provided")
        library_prep_pipeline = LibraryPrepPipeline()
        library_prep_pipeline.execute(configuration, global_configuration, output_directory, log_directory, molecule_packet, rng)

class BeersLibraryPrepValidationException(Exception):
    pass
