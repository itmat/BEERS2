import importlib
import pathlib
import sys
import time
import numpy as np
from beers.cluster_packet import ClusterPacket

class SequencePipeline():
    """
    The class runs all the steps in the sequence pipeline as described and wired together in the configuration
    file. The point of entry into this class is the static method main().
    """

    stage_name = "sequence_pipeline"
    package = "beers.sequence"

    def __init__(self):
        """
        """
        pass

    @staticmethod
    def validate(configuration: dict, global_config: dict):
        """
        Static method to run each step validate process to identify errant parameters.  If any errors are found,
        a BeersSequenceValidationException is raised. Prints error messages to stderr

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
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)

            # Valdiate steps's parameters and print out their errors (if any)
            errors = step_class.validate(parameters, global_config)
            if len(errors) > 0:
                validation_okay = False
                print(f"Validation errors in {step['step_name']}:", file=sys.stderr)
                for error in errors:
                    print('\t' + error, file=sys.stderr)

        if not validation_okay:
            raise BeersSequenceValidationException("Validation error: see stderr for details.")

    def execute(self,
            configuration: dict,
            global_config: dict,
            input_cluster_packet: ClusterPacket,
            output_packet_path: str,
            log_paths: str,
            rng: np.random.Generator,
        ):
        """
        Performs the Sequence Pipeline on a ClusterPacket.

        Parameters
        ----------
        configuration:
            json-like object of the configuration data relevant to this pipeline stage.
        global_config:
            json-like object of the full configuration data
        input_cluster_packet:
            the cluster packet to run through this pipeline stage
        output_packet_path:
            the path to serialize the results to
        log_paths:
            list of paths to write out the logs to, in the same order as each the steps
            appear in the config file
        rng:
            random number generator to use
        """

        # Load and instantiate all steps listed in configuration prior to executing them below.
        steps = []
        for step, log_path in zip(configuration['steps'], log_paths):
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)
            steps.append(step_class(log_path, parameters, global_config))

        print(f"Execution of the {SequencePipeline.stage_name} Started...")
        pipeline_start = time.time()
        cluster_packet = input_cluster_packet
        for step in steps:
            cluster_packet = step.execute(cluster_packet, rng)

        pipeline_elapsed_time = time.time() - pipeline_start
        print(f"Finished sequence pipeline in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a gzip file for inspection
        cluster_packet.serialize(output_packet_path)
        print(f"Output final sample to {output_packet_path}")

    @staticmethod
    def main(
            seed: int,
            configuration: dict,
            global_configuration: dict,
            input_packet_path: str,
            output_packet_path: str,
            log_paths: str
        ):
        """
        Prepares the pipeline, loads the cluster packet, and then executes() the sequence pipeline.

        Parameters
        ----------
        seed:
            value to use as the seed for the random number generator.
        configuration:
            the json-like object containing the configration data specific to the sequencing prep pipeline
        global_configuration:
            json-like object  containing the full configuratino of the entire run
        input_packet_path:
            the file path to the cluster packet to read in
        output_packet_path:
            the file path to output the packet to
        log_paths:
            list of paths to write out the logs to, in the same order as each the steps
            appear in the config file
        """

        cluster_packet = ClusterPacket.deserialize(input_packet_path)
        sample_id = cluster_packet.sample.sample_id
        packet_id = cluster_packet.cluster_packet_id

        # Seed the RNG by the given (global) seed, plus the sample ID, the packet ID, and a constant for the sequence stage of 2
        # This ensures that each packet is repeatable but that each gets their own seed state and that the seed state differs
        # from each of the stages it goes through.
        seed_list = [seed, sample_id, packet_id, 2]
        print(f"Initializing with seed {seed_list}")
        rng = np.random.default_rng(seed_list)

        sequence_pipeline = SequencePipeline()
        sequence_pipeline.execute(configuration, global_configuration, cluster_packet, output_packet_path, log_paths, rng)


class BeersSequenceValidationException(Exception):
    pass
