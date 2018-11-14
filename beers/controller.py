import json
import time
from timeit import default_timer as timer
import numpy as np
from datetime import datetime
from beers.utilities.general_utils import GeneralUtils
from beers.expression.expression_pipeline import ExpressionPipeline
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline

class Controller:

    def run_expression_pipeline(self, args):
        self.retrieve_configuration(args.config)
        self.plant_seed()
        self.create_controller_log()
        ExpressionPipeline.main(self.configuration['expression_pipeline'])

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
        self.seed = self.configuration.get('seed', GeneralUtils.generate_seed())
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