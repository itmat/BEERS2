import subprocess
import os
import json
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from multiprocessing import Pool

class Dispatcher:

    next_molecule_packet_id = 1

    def __init__(self, stage_name, configuration, input_directory_path, output_directory_path):
        self.stage_name = stage_name
        self.configuration = configuration
        self.input_directory_path = input_directory_path
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, "logs")
        self.std_ouput_file_path = os.path.join(log_directory_path, 'std_out.txt')
        self.std_error_file_path = os.path.join(log_directory_path, 'std_err.txt')

    def dispatch(self, method, molecule_packet_filenames):
        if method == 'multicore':
            self.dispatch_multicore(molecule_packet_filenames)
        elif method == 'pmacs':
            self.dispatch_pmacs(molecule_packet_filenames)
        else:
            self.dispatch_serial(molecule_packet_filenames)

    def dispatch_serial(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        stage_process = f"./run_{self.stage_name}.py"
        for molecule_packet_filename in molecule_packet_filenames:
            result = subprocess.call(
            f"{stage_process}"
            f" -c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path}"
            f" -p {molecule_packet_filename}", shell=True)

    def dispatch_multicore(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        data = [(stage_configuration, self.output_directory_path, molecule_packet_filename) for molecule_packet_filename in molecule_packet_filenames]
        pool = Pool(processes=2)
        # TODO - class cannot be hardcoded
        pool.starmap(LibraryPrepPipeline.main, data)

    def dispatch_pmacs(self, molecule_packet_filenames):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        stage_process = f"./run_{self.stage_name}.py"
        for molecule_packet_filename in molecule_packet_filenames:
            # result = subprocess.call(
            #     f"bsub -o {self.std_ouput_file_path} -e {self.std_error_file_path} -J "
            #     f"python ../beers/library_prep/library_prep_pipeline.py"
            #     f" -c '{stage_configuration}' -o {self.output_directory_path}"
            #     f" -p {molecule_packet_filename}", shell=True)
            result = subprocess.call(
                f"bsub -o {self.std_ouput_file_path} -e {self.std_error_file_path} -J "
                f"{stage_process}"
                f" -c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path}"
                f" -p {molecule_packet_filename}", shell=True)