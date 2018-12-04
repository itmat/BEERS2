import importlib
import pickle
import os
import time
import numpy as np
import resource


class SequencePipeline:

    stage_name = "sequence_pipeline"
    package = "beers.sequence"

    def __init__(self, configuration, output_directory_path, cluster_packet):
        self.cluster_packet = cluster_packet
        output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, "logs")
        data_directory_path = os.path.join(output_directory_path, 'data')
        self.log_file_path = os.path.join(log_directory_path, "sequence_pipeline.log")
        self.steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_cluster_pkt{self.cluster_packet.cluster_packet_id}.log"
            step_log_file_path = os.path.join(log_directory_path, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=SequencePipeline.package)
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters))

        results_filename = f"{SequencePipeline.stage_name}_result_cluster_pkt{self.cluster_packet.cluster_packet_id}.pickle"
        self.results_file_path = os.path.join(data_directory_path, results_filename)

    def validate(self, **kwargs):
        if not all([step.validate() for step in self.steps]):
            raise BeersSequenceValidationException("Validation error in step: see stderr for details.")

    def execute(self):
        print(f"Execution of the {SequencePipeline.stage_name} Started...")
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

        # Write final sample to a gzip file for inspection
        cluster_packet.serialize(self.results_file_path)
        print(f"Output final sample to {self.results_file_path}")

    @staticmethod
    def main(configuration, output_directory_path, cluster_packet):
        sequence_pipeline = SequencePipeline(configuration, output_directory_path, cluster_packet)
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass
