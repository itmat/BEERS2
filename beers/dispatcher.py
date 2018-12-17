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
    """
    This is meant to be a singleton object instantiated by the controller.  It creates and executes the subprocess
    commands needed to spawn the library_prep and sequence pipeline processes for each packet (molecule or cluster).
    """

    packet_pattern = re.compile(r'^.*_pkt(\d+).gzip$')

    def __init__(self,
                 run_id,
                 dispatcher_mode,
                 seed,
                 stage_name,
                 configuration,
                 input_directory_path,
                 output_directory_path,
                 directory_structure):
        """
        Because the processes have no knowledge of the controller's state, a lot of information must be conveyed on the
        command line when calling each pipeline stage process.  Consequently the dispatcher object takes a lot of
        parameters.
        :param run_id: The run id the user has designated for this simulation
        :param dispatcher_mode: Three modes are currently supported:  serial, multicore, and lsf
        :param seed: To achieve reproducibility, a seed either provided by the user or the controller is passed along
        to each process.
        :param stage_name: The pipeline stage for which a process is called (currently library_prep_pipeline or
        sequence_pipeline.  The stage_name is used also to build the name of the executable process that resides in
        the bin directory.
        :param configuration: the configuration json data from the configuration file.
        :param input_directory_path: The directory containing the molecule/cluster packets to feed into the pipeline.
        :param output_directory_path: The directory containing the output data from the pipeline
        :param directory_structure: Because of the number of files that can possibly be generated, files are
        partitioned across sub-directories and in the case of enough files, possibly across sub-sub directories.  This
        string representation identifies the number of sub-directory levels and the number of directories on each level.
        This information must be forwarded to each called process.
        """
        # TODO input directory path may now be redundant since the full packet path is sent to the dispatch method
        self.run_id = run_id
        self.dispatcher_mode = dispatcher_mode
        self.seed = seed
        self.stage_name = stage_name
        self.configuration = configuration
        self.input_directory_path = input_directory_path
        self.output_directory_path = output_directory_path
        self.log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        self.directory_structure = directory_structure

    def dispatch(self, packet_file_paths):
        """
        This is the entry method for dispatching jobs (calls to library prep or sequence pipelines).  A different
        dispatch method is called depending on the dispatcher mode:  serial, lsf, or multicore.
        :param packet_file_paths: The absolute file paths to the input packets
        """
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
            command = f"{stage_process} " \
                      f"-s {self.seed} " \
                      f"-c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path} " \
                      f"-p {packet_file_path} -d {self.directory_structure}"
            subprocess.call(command, shell=True)

    def dispatch_multicore(self, packet_file_paths):
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        data = [(self.seed, stage_configuration, self.output_directory_path, packet_file_path, self.directory_structure)
                for packet_file_path in packet_file_paths]
        pool = Pool(processes=2)
        if self.stage_name == 'library_prep_pipeline':
            pool.starmap(LibraryPrepPipeline.main, data)
        else:
            pool.starmap(SequencePipeline.main, data)

    def dispatch_lsf(self, packet_file_paths):
        """
        The method dispatches to various nodes, one job for each packet as indicated by the packet file path list.  The
        job name is a combination of the run id, the stage name and the packet id, which should hopefully insure
        uniqueness.
        :param packet_file_paths:  The absolute file paths to the input packets
        """
        stage_configuration = json.dumps(self.configuration[self.stage_name])
        stage_process = f"./run_{self.stage_name}.py"
        for packet_file_path in packet_file_paths:
            stdout_file_path, stderr_file_path, packet_id = self.get_stdout_and_stderr_subdirectories(packet_file_path)
            command = f"bsub -o {stdout_file_path} -e {stderr_file_path} " \
                      f"-J run{self.run_id}_{self.stage_name}_pkt{packet_id} " \
                      f"{stage_process} " \
                      f"-s {self.seed} " \
                      f"-c '{stage_configuration}' -i {self.input_directory_path} -o {self.output_directory_path} " \
                      f"-p {packet_file_path} -d {self.directory_structure}"
            subprocess.call(command, shell=True)

    @staticmethod
    def get_packet_id_from_file(packet_file_path):
        """
        Helper method to retrieve the packet id from the packet's file path.  The packet filename is expected to end
        with _pkt\d+.gzip, which we use here to extract the packet name.  Since the software creates all packets, it is
        expected that any failure here is a software error.
        :param packet_file_path: The packet file path from which the packet id is extracted.
        :return: The packet id cast into an integer
        """
        packet_filename = os.path.split(packet_file_path)[1]
        packet_pattern_match = re.match(Dispatcher.packet_pattern, packet_filename)
        if packet_pattern_match:
            return int(packet_pattern_match.group(1))
        assert False, f"The packet filename, {packet_filename} does not follow convention."

    def get_stdout_and_stderr_subdirectories(self, packet_file_path):
        """
        For LSF, stdout/stderr filenames must be created such that there is no collision between those generated by
        other jobs.  Since each job handles one packet, the packet's packet id is affixed to the end of each stdout/
        stderr file combination.  Additionally, these outputs are organized into sub-directories in the same manner
        as are the various log and data outputs of the pipelines because here too, the number of files generated may
        exceed the number that can reasonably kept in one directory.
        :param packet_file_path: The packet file path from which the packet id is extracted.
        :return: a tuple containing the full paths of the stdout and stderr, in that order
        """
        packet_id = Dispatcher.get_packet_id_from_file(packet_file_path)
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
