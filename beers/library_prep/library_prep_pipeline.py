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

class LibraryPrepPipeline():
    """
    The class runs all the steps in the library_prep pipeline as described and wired together in the configuration
    file.  The point of entry into this class is the static method main().
    """

    stage_name = "library_prep_pipeline"
    pipeline_log_subdirectory_name = "Pipeline"
    package = "beers.library_prep"

    def __init__(self):
        pass

    @staticmethod
    def validate(configuration: dict, global_config: dict):
        """
        Static method to run each step validate process to identify errant parameters.  If any errors are found,
        a BeersLibraryPrepValidationException is raised. Prints error messages to stderr.

        Parameters
        ----------
        configuration:
            json-like object containing the SequencePipeline-specific configuration
        global_config:
            json-like object containing the full config of the BEERS run
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

    def execute(self,
            configuration: dict,
            global_config: dict,
            output_directory: str,
            log_directory: str,
            input_molecule_packet: MoleculePacket,
            rng: np.random.Generator,
        ):
        """
        Performs the Sequence Pipeline on a MoleculePacket.

        Parameters
        ----------
        configuration:
            json-like object of the configuration data relevant to this pipeline stage.
        global_config:
            json-like object of the full configuration data
        output_directory:
            the folder to serialize the results to as 'library_prep_pipeline_results_molecule_pkt###.txt'
            as well as the .quant_file for both results and inputs.
        log_directory:
            directory to output logs to
        input_molecule_packet:
            the molecule packet to run through this pipeline stage
        rng:
            random number generator to use
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

    def print_summary(self, original_ids: set[int], sample: list[Molecule], elapsed_time:float=None):
        """Output a summary of the sample (number of molecules, time taken, etc.).

        original_ids:
            Molecule IDs from input molecule packet
        sample:
            list of molecules in the sample
        elapsed_time:
            elapsed time to display, if not None
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
    def main(
            seed: int,
            configuration: dict,
            global_configuration: dict,
            output_directory: str,
            log_directory: str,
            molecule_packet_filename: str,
            packet_id: str,
            distribution_directory: str = None,
            molecules_per_packet_from_distribution: int = 10000,
            sample_id: str= None,
            ):
        """
        This method would be called by a command line script in the bin directory.  It sets a random seed, loads a
        directory containing the relevant parts of the user's configuration file, unmarshalls a molecule packet from
        the provided molecule packet filename, initializes and validates the library prep pipeline stage and then
        executes it for the molecule packet.

        Parameters
        ----------
        seed:
            value to use as the seed for the random number generator
        configuration:
            the json string containing the configration data specific to the library prep pipeline
        global_configuration:
            json string containing the configuration data for the whole run
        output_directory:
            path to directory where we will output molecule files and quant files
        log_directory:
            path to log directories
        directory_structure:
            instructions for creating the scaffolding needed to house the pipeline data and logs
        molecule_packet_filename:
            the file from which to unmarshall the molecule packet
        packet_id:
            id number to assign the packet
        distribution_directory:
            directory of CAMPAREE output distribution datas from which to generate molecules
            (default None, must supply molecule_packet_filename instead)
        molecules_per_packet_from_distribution:
            packet size to generate if using distributions
        sample_id:
            id number of the sample (if None, derive from the molecule packet provided)
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
