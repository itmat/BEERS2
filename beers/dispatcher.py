import subprocess
import os
import re
import json
from beers.utilities.general_utils import GeneralUtils
from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline
from beers.sequence.sequence_pipeline import SequencePipeline
from multiprocessing import Pool
from beers.constants import CONSTANTS


class Dispatcher:

    packet_pattern = '^.*_pkt(\d+).gzip$'

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
        self.log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
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
        for packet_file_path in packet_file_paths:
            stdout_file_path, stderr_file_path, packet_id = self.get_stdout_and_stderr_subdirectories(packet_file_path)
            command = f"bsub -o {stdout_file_path} -e {stderr_file_path} -J {self.stage_name}_pkt{packet_id} " \
                f"{stage_process} " \
                f"-c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path} " \
                f"-p {packet_file_path} -d {self.directory_structure}"
            print(command)
            result = subprocess.call(command, shell=True)

    def get_packet_id_from_file(self, packet_file_path):
        packet_filename = os.path.split(packet_file_path)[1]
        packet_pattern_match = re.match(Dispatcher.packet_pattern, packet_filename)
        if packet_pattern_match:
            return int(packet_pattern_match.group(1))
        assert False, f"The packet filename, {packet_filename} does not follow convention."

    def get_stdout_and_stderr_subdirectories(self, packet_file_path):
        packet_id = self.get_packet_id_from_file(packet_file_path)
        subdirectory_list = GeneralUtils.get_output_subdirectories(packet_id, self.directory_structure)
        stdout_filename = f"{CONSTANTS.STDOUT_SUBDIRECTORY_NAME}_pkt{packet_id}.log"
        stdout_file_path = os.path.join(self.log_directory_path,
                                        CONSTANTS.STDOUT_SUBDIRECTORY_NAME,
                                        *subdirectory_list,
                                        stdout_filename)
        stderr_filename = f"{CONSTANTS.STDERR_SUBDIRECTORY_NAME}_pkt{packet_id}.log"
        stderr_file_path = os.path.join(self.log_directory_path,
                                        CONSTANTS.STDERR_SUBDIRECTORY_NAME,
                                        *subdirectory_list,
                                        stderr_filename)
        return stdout_file_path, stderr_file_path, packet_id