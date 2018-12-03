import subprocess
import os
import json
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from multiprocessing import Pool

class Dispatcher:

    next_molecule_packet_id = 1

    def __init__(self, stage_name, configuration, output_directory_path):
        self.stage_name = stage_name
        self.configuration = configuration
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, "logs")
        self.std_ouput_file_path = os.path.join(log_directory_path, 'std_out.txt')
        self.std_error_file_path = os.path.join(log_directory_path, 'std_err.txt')

    def dispatch(self, method, molecule_packet_filenames):
        if method == 'multicore':
            self.dispatch_multicore(molecule_packet_filenames)
        elif method == 'pmcas':
            self.dispatch_pmacs(molecule_packet_filenames)
        else:
            self.dispatch_serial(molecule_packet_filenames)

    def dispatch_serial(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration['library_prep_pipeline'])
        for molecule_packet_filename in molecule_packet_filenames:
            result = subprocess.call(
            f"python ../beers/library_prep/library_prep_pipeline.py"
            f" -c '{stage_configuration}' -o {self.output_directory_path}"
            f" -p {molecule_packet_filename}", shell=True)

    def dispatch_multicore(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration['library_prep_pipeline'])
        data = [(stage_configuration, self.output_directory_path, molecule_packet_filename) for molecule_packet_filename in molecule_packet_filenames]
        pool = Pool(processes=2)
        pool.starmap(LibraryPrepPipeline.main, data)

    def dispatch_pmacs(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration['library_prep_pipeline'])
        for molecule_packet_filename in molecule_packet_filenames:
            result = subprocess.call(
                f"bsub -o {self.std_ouput_file_path} -e {self.std_error_file_path} -J "
                f"python ../beers/library_prep/library_prep_pipeline.py"
                f" -c '{stage_configuration}' -o {self.output_directory_path}"
                f" -p {molecule_packet_filename}", shell=True)
