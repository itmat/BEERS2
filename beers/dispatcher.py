import subprocess
import os
import re
import json
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline
from multiprocessing import Pool

class Dispatcher:

    next_molecule_packet_id = 1

    def __init__(self,
                 dispatcher_mode,
                 stage_name,
                 configuration,
                 input_directory_path,
                 output_directory_path,
                directory_structure):
        self.dispatcher_mode = dispatcher_mode
        self.stage_name = stage_name
        self.configuration = configuration
        self.input_directory_path = input_directory_path
        self.output_directory_path = output_directory_path
        self.log_directory_path = os.path.join(output_directory_path, "logs")
        self.directory_structure = directory_structure

    def dispatch(self, packet_file_paths):
        if self.dispatcher_mode == 'multicore':
            self.dispatch_multicore(packet_file_paths)
        elif self.dispatcher_mode == 'lsf':
            self.dispatch_lsf(packet_file_paths)
        else:
            self.dispatch_serial(packet_file_paths)

    def dispatch_serial(self, packet_file_paths):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        stage_process = f"./run_{self.stage_name}.py"
        for packet_file_path in packet_file_paths:
            result = subprocess.call(
            f"{stage_process}"
            f" -c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path}"
            f" -p {packet_file_path} -d {self.directory_structure}", shell=True)

    def dispatch_multicore(self, packet_file_paths):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        data = [(stage_configuration, self.output_directory_path, packet_file_path, self.directory_structure)
                for packet_file_path in packet_file_paths]
        pool = Pool(processes=2)
        if self.stage_name == 'library_prep_pipeline':
            pool.starmap(LibraryPrepPipeline.main, data)
        else:
            pool.starmap(SequencePipeline.main, data)

    def dispatch_lsf(self, packet_file_paths):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        stage_process = f"./run_{self.stage_name}.py"
        for ctr, packet_file_path in enumerate(packet_file_paths):
            std_ouput_file_path = os.path.join(self.log_directory_path, f'std_out_{ctr}.txt')
            std_error_file_path = os.path.join(self.log_directory_path, f'std_err_{ctr}.txt')
            result = subprocess.call(
                f"bsub -o {std_ouput_file_path} -e {std_error_file_path} -J {self.stage_name}_{ctr} "
                f"{stage_process}"
                f" -c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path}"
                f" -p {packet_file_path} -d {self.directory_structure}", shell=True)