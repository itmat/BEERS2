import json
import time
from timeit import default_timer as timer
import numpy as np
import os
from datetime import datetime
from beers.utilities.general_utils import GeneralUtils
from beers.expression.expression_pipeline import ExpressionPipeline
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline
from beers.sample import Sample

class Controller:

    def run_expression_pipeline(self, args):
        self.retrieve_configuration(args.config)
        self.plant_seed()
        self.create_controller_log()
        self.assemble_input_samples()
        ExpressionPipeline.main(self.input_samples, self.configuration['expression_pipeline'])

    def run_library_prep_pipeline(self, args):
        self.retrieve_configuration(args.config)
        self.plant_seed()
        self.create_controller_log()
        LibraryPrepPipeline.main(self.configuration['library_prep_pipeline'])

    def run_sequence_pipeline(self, args):
        start = timer()
        self.retrieve_configuration(args.config)
        self.plant_seed()
        self.create_controller_log()
        SequencePipeline.main(self.configuration['sequence_pipeline'])
        end = timer()
        print(f"Sequence Pipeline: {end - start}")

    def retrieve_configuration(self, configuration_file_path):
        with open(configuration_file_path, "r+") as configuration_file:
            self.configuration = json.load(configuration_file)

    def plant_seed(self):
        self.seed = self.configuration['controller'].get('seed', GeneralUtils.generate_seed())
        np.random.seed(self.seed)

    def create_controller_log(self):
        log_file_path = self.configuration['controller']['log_file_path']
        with open(log_file_path, 'w') as controller_log_file:
            timestamp = time.time()
            current_datetime = datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')
            controller_log_file.write(f"Run Timestamp (UTC):\t{current_datetime}\n")
            controller_log_file.write(f"Seed:\t{self.seed}\n")
            controller_log_file.write(f"Configuration:\n")
            json.dump(self.configuration, controller_log_file, indent=2)

    def assemble_input_samples(self):
        input_directory_path = self.configuration['controller']["input"]["directory_path"]
        self.input_samples = []
        for input_sample in self.configuration['controller']["input"]["data"]:
            sample_name = os.path.splitext(input_sample["filename"])[0]
            input_sample_file_path = os.path.join(input_directory_path, input_sample["filename"])
            self.input_samples.append(
                Sample(Sample.next_sample_id,
                       sample_name,
                       input_sample_file_path,
                       input_sample.get("gender", None)))
            Sample.next_sample_id += 1