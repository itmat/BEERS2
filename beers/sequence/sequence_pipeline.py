import argparse
import importlib
import pathlib
import os
import time
import numpy as np
import resource
import json
import re
from beers_utils.constants import CONSTANTS
from beers.cluster_packet import ClusterPacket
from beers_utils.general_utils import GeneralUtils

from beers.abstract_beers_pipeline import AbstractBeersPipeline

class SequencePipeline():
    """
    The class runs all the steps in the sequence pipeline as described and wired together in the configuration
    file.  The point of entry into this class is the static method main().
    """

    stage_name = "sequence_pipeline"
    package = "beers.sequence"

    def __init__(self):
        """
        Many initialization steps here include identifying the log and data directories and subdirectories so
        that data and log files generated are placed in the correct locations.  Note that the directory structure
        has already been created by the controller.  The steps described in the configuration dictionary are
        instantiated and those instantialed steps are added to a list to further use.
        """
        pass

    @staticmethod
    def validate(configuration, global_configuration):
        """
        Static method to run each step validate process to identify errant parameters.  If any errors are found,
        a validation exception is raised.
        """
        steps = []
        # Find the step classes in the configuration
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)
            steps.append(step_class(None, parameters, global_configuration))
        # Validate configuration of each step
        if not all([step.validate() for step in steps]):
            raise BeersSequenceValidationException("Validation error in step: see stderr for details.")

    def execute(self, configuration, global_config, input_cluster_packet, output_packet_path, log_directory, rng):
        """
        Opens the pipeline log for writing and serially runs the execute method of each step object found in the
        step list generated when this pipeline stage was initialized.  The final product (a cluster packet
        modified with additional information) is serialized into a data file.
        :param configuration: dictionary of the configuration data relevant to this pipeline stage.
        :param global_config:  dictionary of the full configuration data
        :param input_cluster_packet: the cluster packet to run through this pipeline stage
        :param output_packet_path: the path to serialize the results to
        :param log_directory: directory to output logs to
        :param rng: rnadom number generator to use
        :return:
        """

        log_directory = pathlib.Path(log_directory)
        log_file_path = log_directory / f"sequence_pipeline_cluster_pkt{input_cluster_packet.cluster_packet_id}.log"

        # Load and instantiate all steps listed in configuration prior to executing them below.
        self.steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_cluster_pkt{input_cluster_packet.cluster_packet_id}.log"
            step_log_dir = log_directory / step_name
            step_log_dir.mkdir(exist_ok=True)
            step_log_file_path = step_log_dir / step_log_filename
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters, global_config))

        print(f"Execution of the {SequencePipeline.stage_name} Started...")
        with open(log_file_path, 'w') as log_file:
            pipeline_start = time.time()
            cluster_packet = input_cluster_packet
            for step in self.steps:
                cluster_packet = step.execute(cluster_packet, rng)

            pipeline_elapsed_time = time.time() - pipeline_start
            print(f"Finished sequence pipeline in {pipeline_elapsed_time:.1f} seconds")

            # Write final sample to a gzip file for inspection
            cluster_packet.serialize(output_packet_path)
            print(f"Output final sample to {output_packet_path}")
            log_file.write("Sequencing pipeline completed successfully\n")

    @staticmethod
    def main(seed, configuration, global_configuration, input_packet_path, output_packet_path, log_directory):
        """
        This method would be called by a command line script in the bin directory.  It sets a random seed, loads a
        directory containing the relevant parts of the user's configuration file, unmarshalls a cluster packet from
        the provided cluster packet filename, initializes and validates the sequence pipeline stage and then
        executes it for the cluster packet.  Since the controller cannot conclude until all sequence pipelines are
        run, the last action taken by the sequence pipeline is to note its completion to the auditor, which records
        the cluster id in an audit file.  This happens regardless of the outcome of this pipeline stage.
        :param seed: value to use as the seed for the random number generator
        :param configuration: the json string containing the configration data specific to the sequencing prep pipeline
        :param global_configuration: json string containing the full configuratino of the entire run
        :param input_packet_path: the file path to the cluster packet to read in
        :param output_packet_path: the file path to output the packet to
        :param log_directory: directory path to output logs to
        """

        # Normally the cluster_packet_id should be derived from the serialized cluster_packet in the file.  But if
        # the file cannot be found, we still need an id to report back to the auditor if at all possible.  So we
        # extract it from the file name just in case.
        cluster_packet = None
        cluster_packet_filename = str(pathlib.Path(input_packet_path).name)
        cluster_packet_id_pattern = re.compile(r'^.*cluster_packet.*_pkt(\d+)\..*$')

        cluster_packet = ClusterPacket.get_serialized_cluster_packet(input_packet_path)
        sample_id = cluster_packet.sample.sample_id
        packet_id = cluster_packet.cluster_packet_id

        # Seed the RNG by the given (global) seed, plus the sample ID, the packet ID, and a constant for the sequence stage of 2
        # This ensures that each packet is repeatable but that each gets their own seed state and that the seed state differs
        # from each of the stages it goes through.
        seed_list = [seed, sample_id, packet_id, 2]
        print(f"Initializing with seed {seed_list}")
        rng = np.random.default_rng(seed_list)

        sequence_pipeline = SequencePipeline()
        sequence_pipeline.execute(configuration, global_configuration, cluster_packet, output_packet_path, log_directory, rng)


class BeersSequenceValidationException(Exception):
    pass
