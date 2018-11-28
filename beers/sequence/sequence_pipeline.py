import importlib
import pickle
import os
import time
import numpy as np
import resource


class SequencePipeline:

    def __init__(self, configuration, cluster_packet=None):
        self.steps = []
        input_directory_path = configuration["input"]["directory_path"]
        output_directory_path = configuration["output"]["directory_path"]
        self.log_file_path = os.path.join(output_directory_path, configuration["output"]["log_filename"])
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = step["log_filename"]
            step_log_file_path = os.path.join(output_directory_path, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package="beers.sequence")
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters))
        self.cluster_packet = cluster_packet
        if not self.cluster_packet:
            cluster_packet_filename = configuration["input"]["cluster_packet_filename"]
            self.molecule_packet_file_path = os.path.join(input_directory_path, cluster_packet_filename)
            self.populate_cluster_packet()
        results_filename = configuration["output"]["results_filename"]
        self.results_file_path = os.path.join(output_directory_path, results_filename)

    def validate(self, **kwargs):
        if not all([step.validate() for step in self.steps]):
            raise BeersSequenceValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        print("Execution of the Sequence Pipeline Started...")
        with open(self.log_file_path, 'w') as log_file:
            pipeline_start = time.time()
            cluster_packet = self.cluster_packet
            for step in self.steps:
                step_start = time.time()
                cluster_packet = step.execute(cluster_packet)
                elapsed_time = time.time() - step_start

                random_state = np.random.get_state()
                log_file.write(f"# random state following {step.__class__.__name__} is {random_state}\n")
                print(f"{step.__class__.name} complete - process RAM currently at"
                      f" {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")

            pipeline_elapsed_time = time.time() - pipeline_start
            print(f"Finished sequence pipeline in {pipeline_elapsed_time:.1f} seconds")

        # Write final sample to a pickle file for inspection
        with open(self.results_file_path, "wb") as results_file:
            pickle.dump(cluster_packet, results_file)
        print(f"Output final sample to {self.results_file_path}")

    def populate_cluster_packet(self):
        print(f"Reading molecule packet from pickle file {self.cluster_packet_file_path}")
        with open(self.cluster_packet_file_path, 'rb') as cluster_packet_file:
            self.cluster_packet = pickle.load(cluster_packet_file)

    @staticmethod
    def main(configuration, cluster_packet=None):
        sequence_pipeline = SequencePipeline(configuration, cluster_packet)
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass

