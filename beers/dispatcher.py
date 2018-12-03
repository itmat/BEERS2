import subprocess
import os
import json
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from multiprocessing import Pool

class Dispatcher:

    next_molecule_packet_id = 1

    def __init__(self, configuration):
        self.configuration = configuration

    def dispatch(self, method, output_directory_path, molecule_packet_filenames):
        if method == 'multicore':
            self.dispatch_multicore(output_directory_path, molecule_packet_filenames)
        elif method == 'pmcas':
            self.dispatch_pmacs(output_directory_path, molecule_packet_filenames)
        else:
            self.dispatch_serial(output_directory_path, molecule_packet_filenames)

    def dispatch_serial(self, output_directory_path, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration['library_prep_pipeline'])
        for molecule_packet_filename in molecule_packet_filenames:
            result = subprocess.call(
            f"python ../beers/library_prep/library_prep_pipeline.py"
            f" -c '{stage_configuration}' -o {output_directory_path}"
            f" -p {molecule_packet_filename}", shell=True)

    def dispatch_multicore(self, output_directory_path, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration['library_prep_pipeline'])
        data = [(stage_configuration, output_directory_path, molecule_packet_filename) for molecule_packet_filename in molecule_packet_filenames]
        pool = Pool(processes=2)
        pool.starmap(LibraryPrepPipeline.main, data)

    def dispatch_pmacs(self, output_directory_path, molecule_packet_filenames):
        pass
