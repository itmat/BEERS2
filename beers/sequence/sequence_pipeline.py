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
            log_directory: str,
            rng: np.random.Generator,
        ):
        """
        Performs the Sequence Pipeline on a ClusterPacket.

        Parameters
        ----------
        configuration:
            dictionary of the configuration data relevant to this pipeline stage.
        global_config:
            dictionary of the full configuration data
        input_cluster_packet:
            the cluster packet to run through this pipeline stage
        output_packet_path:
            the path to serialize the results to
        log_directory:
            directory to output logs to
        rng:
            random number generator to use
        """

        log_directory = pathlib.Path(log_directory)
        log_file_path = log_directory / f"sequence_pipeline_cluster_pkt{input_cluster_packet.cluster_packet_id}.log"

        # Load and instantiate all steps listed in configuration prior to executing them below.
        steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_cluster_pkt{input_cluster_packet.cluster_packet_id}.log"
            step_log_dir = log_directory / step_name
            step_log_dir.mkdir(exist_ok=True)
            step_log_file_path = step_log_dir / step_log_filename
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)
            steps.append(step_class(step_log_file_path, parameters, global_config))

        print(f"Execution of the {SequencePipeline.stage_name} Started...")
        with open(log_file_path, 'w') as log_file:
            pipeline_start = time.time()
            cluster_packet = input_cluster_packet
            for step in steps:
                cluster_packet = step.execute(cluster_packet, rng)

            pipeline_elapsed_time = time.time() - pipeline_start
            print(f"Finished sequence pipeline in {pipeline_elapsed_time:.1f} seconds")

            # Write final sample to a gzip file for inspection
            cluster_packet.serialize(output_packet_path)
            print(f"Output final sample to {output_packet_path}")
            log_file.write("Sequencing pipeline completed successfully\n")

    @staticmethod
    def main(
            seed: int,
            configuration: dict,
            global_configuration: dict,
            input_packet_path: str,
            output_packet_path: str,
            log_directory: str
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
        log_directory:
            directory path to output logs to
        """

        # Normally the cluster_packet_id should be derived from the serialized cluster_packet in the file.  But if
        # the file cannot be found, we still need an id to report back to the auditor if at all possible.  So we
        # extract it from the file name just in case.
        cluster_packet = None

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
        sequence_pipeline.execute(configuration, global_configuration, cluster_packet, output_packet_path, log_directory, rng)


class BeersSequenceValidationException(Exception):
    pass
