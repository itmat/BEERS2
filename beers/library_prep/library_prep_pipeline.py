import argparse
import importlib
import pathlib
import time
import re
import os
import resource
import json
import math

import numpy as np

from beers_utils.sample import Sample
from beers_utils.constants import CONSTANTS
from beers_utils.molecule import Molecule
from beers_utils.molecule_packet import MoleculePacket
from beers_utils.general_utils import GeneralUtils

from camparee.molecule_maker import MoleculeMakerStep

from beers.abstract_beers_pipeline import AbstractBeersPipeline

class LibraryPrepPipeline(AbstractBeersPipeline):
    """
    The class runs all the steps in the library_prep pipeline as described and wired together in the configuration
    file.  The point of entry into this class is the static method main().  In HPC scenarios, this object's main is
    invoked on each node included in the compute, with a different molecule packet delivered to each node.
    """

    stage_name = "library_prep_pipeline"
    pipeline_log_subdirectory_name = "Pipeline"
    package = "beers.library_prep"

    def __init__(self, output_directory_path):
        """
        Many initialization steps here include identifying the log and data directories and subdirectories so
        that data and log files generated are placed in the correct locations.
        :param output_directory_path: top level directory for data and logs generated by this pipeline stage
        """
        self.log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        self.data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)

        """
        Moving molecule handling to execute method. This means BEERS2 will instantiate a
        single LibraryPrepPipeline object which contains all of the configurations common
        to all potential pipelines (executed on different molecule packets). This is to make
        BEERS2 conform to the same job handling protocol as BEERS_UTILS and CAMPAREE. It
        does break with the OO convention of making an instance of the library prep pipeline
        for each molecule packet. This change may be temporary, and we determine it's better
        to update BEERS_UTILS and CAMPAREE to more closely follow an object-oriented design.

        Keeping the original code here, for now, so it's easier to revert the change if we
        decide to do so. We can remove this comment once we're happy with the implementation.

        self.molecule_packet = molecule_packet
        self.subdirectory_list = \
            GeneralUtils.get_output_subdirectories(self.molecule_packet.molecule_packet_id, directory_structure)
        data_subdirectory_path = os.path.join(data_directory_path, *subdirectory_list)
        self.original_ids = set(str(m.molecule_id) for m in self.molecule_packet.molecules)
        self.print_summary(self.molecule_packet.molecules)
        self.log_file_path = os.path.join(self.log_directory_path,
                                          LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                          *subdirectory_list,
                                          f"{LibraryPrepPipeline.stage_name}_"
                                          f"molecule_pkt{self.molecule_packet.molecule_packet_id}.log")
        self.quant_file_path = os.path.join(self.log_directory_path,
                                          LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                          *subdirectory_list,
                                          f"{LibraryPrepPipeline.stage_name}_"
                                          f"molecule_pkt{self.molecule_packet.molecule_packet_id}.quant_file")
        self.global_config = global_config
        self.steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_molecule_pkt{self.molecule_packet.molecule_packet_id}.log"
            step_log_file_path = os.path.join(self.log_directory_path, step_name, *subdirectory_list, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)
            self.steps.append(step_class(step_log_file_path, parameters, global_config))
        results_filename = f"{LibraryPrepPipeline.stage_name}_" \
                           f"result_molecule_pkt{self.molecule_packet.molecule_packet_id}.txt"
        self.results_file_path = os.path.join(data_subdirectory_path, results_filename)
        results_quant_filename = f"{LibraryPrepPipeline.stage_name}_" \
                           f"result_molecule_pkt{self.molecule_packet.molecule_packet_id}.quant_file"
        self.results_quant_file_path = os.path.join(data_subdirectory_path, results_quant_filename)
        """

    @staticmethod
    def validate(configuration, global_config):
        """
        Static method to run each step validate process to identify errant parameters.  If any errors are found,
        a validation exception is raised.
        """
        #if not all([molecule.validate() for molecule in self.molecule_packet.molecules]):
        #    raise BeersLibraryPrepValidationException("Validation error in molecule packet: see stderr for details.")
        steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)
            steps.append(step_class(None, parameters, global_config))
        if not all([step.validate() for step in steps]):
            raise BeersLibraryPrepValidationException("Validation error in step: see stderr for details.")

    def execute(self, configuration, global_config, directory_structure, input_molecule_packet):
        """
        Opens the pipeline log for writing and serially runs the execute method of each step object found in the
        step list generated when this pipeline stage was initialized.  The final product (a modified molecule
        packet) is serialized into a data file.
        Note that the directory structure
        has already been created by the controller.  The steps described in the configuration dictionary are
        instantiated and those instantiated steps are added to a list for later invocation.
        :param configuration:  dictionary of the configuration data relevant to this pipeline stage.
        :param global_config: dictionary of full config data
        :param directory_structure: instructions for navigating the subdirectories under the output directory
        :param input_molecule_packet: the molecule packet to run through this pipeline stage
        """

        subdirectory_list = \
            GeneralUtils.get_output_subdirectories(input_molecule_packet.molecule_packet_id, directory_structure)
        data_subdirectory_path = os.path.join(self.data_directory_path, *subdirectory_list)
        original_ids = set(str(m.molecule_id) for m in input_molecule_packet.molecules)
        self.print_summary(original_ids, input_molecule_packet.molecules)
        log_file_path = os.path.join(self.log_directory_path,
                                     LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                     *subdirectory_list,
                                     f"{LibraryPrepPipeline.stage_name}_"
                                     f"molecule_pkt{input_molecule_packet.molecule_packet_id}.log")
        quant_file_path = os.path.join(self.log_directory_path,
                                       LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                       *subdirectory_list,
                                       f"{LibraryPrepPipeline.stage_name}_"
                                       f"molecule_pkt{input_molecule_packet.molecule_packet_id}.quant_file")

        results_filename = f"{LibraryPrepPipeline.stage_name}_" \
                           f"result_molecule_pkt{input_molecule_packet.molecule_packet_id}.txt"
        results_file_path = os.path.join(data_subdirectory_path, results_filename)
        results_quant_filename = f"{LibraryPrepPipeline.stage_name}_" \
                                 f"result_molecule_pkt{input_molecule_packet.molecule_packet_id}.quant_file"
        results_quant_file_path = os.path.join(data_subdirectory_path, results_quant_filename)

        # Initialize steps prior to executing them below.
        pipeline_steps = []
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}_molecule_pkt{input_molecule_packet.molecule_packet_id}.log"
            step_log_file_path = os.path.join(self.log_directory_path, step_name, *subdirectory_list, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package=LibraryPrepPipeline.package)
            step_class = getattr(module, step_name)
            pipeline_steps.append(step_class(step_log_file_path, parameters, global_config))

        # TODO what about gzipping?
        with open(log_file_path, 'w') as log_file:
            # self.log_sample(log_file)
            # Log the molecules the pipeline stage begins with
            log_file.write(Molecule.header)
            for molecule in input_molecule_packet.molecules:
                log_file.write(molecule.log_entry())
            pipeline_start = time.time()
            molecule_packet = input_molecule_packet
            molecule_packet.write_quantification_file(quant_file_path)
            for step in pipeline_steps:
                step_name = step.__class__.name if hasattr(step.__class__, 'name') else step.__class__.__name__
                step_start = time.time()
                molecule_packet = step.execute(molecule_packet)
                elapsed_time = time.time() - step_start

                self.print_summary(original_ids, molecule_packet.molecules, elapsed_time)

                # Saves the state of the random number generator after completion of each step in the event that
                # we may need to do to partial repeat of the library_prep pipeline stage (although exactly how
                # to implement such a partial repeat hasn't been given serious consideration as yet).
                random_state = np.random.get_state()
                log_file.write(f"# random state following {step_name} is {random_state}\n")
                print(f"{step_name} complete - process RAM currently at"
                      f" {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6} GB")

            pipeline_elapsed_time = time.time() - pipeline_start
            print(f"Finished {LibraryPrepPipeline.stage_name} in {pipeline_elapsed_time:.1f} seconds")

            # Write final sample to a gzip file for inspection
            molecule_packet.serialize(results_file_path)
            print(f"Output final sample to {results_file_path}")
            molecule_packet.write_quantification_file(results_quant_file_path)
            print(f"Output final sample quantification to {results_quant_file_path}")
            log_file.write("Library prep pipeline completed successfully\n")

    def get_commandline_call(self, configuration, global_config, directory_structure, input_molecule_packet):
        """
        Prepare command to execute the LibraryPrepPipeline from the command line,
        given all of the arguments used to run the execute() function.

        NOTE: This function currently does nothing and is included to satisfy the requirements
              of the AbstractBeersPipeline interface, which this class implements. At the
              moment, the command line call for this pipeline is assembled by the dispatcher,
              without any input from the class itself. This code is an intermediate step along
              the way to updating BEERS2 to fully support the BEERS_UTILS job submission
              framework. At the moment, the main function transforms a bunch of the command
              line arguments into objects, before passing them to the execute function. This
              makes it impossible to retrieve the original text of the command line arguments
              from these objects.

        """

        # Retrieve path to the library_prep_pipeline.py script.
        library_prep_pipeline_path = os.path.realpath(__file__)
        # If the above command returns a string with a "pyc" extension, instead
        # of "py", strip off "c" so it points to this script.
        library_prep_pipeline_path = library_prep_pipeline_path.rstrip('c')

        command = (f" python {library_prep_pipeline_path}")
        return command

    """
    Current form of this class does not have a self.molecule_package object.
    Since log_sample method was only called once, moved this code to the
    execute method.

    def log_sample(self, log_file):
        \"""
        Logs the molecules the pipeline stage begins with.  We need to
        consider gzipping.
        :param log_file: The handle to the log file.
        \"""
        log_file.write(Molecule.header)
        for molecule in self.molecule_packet.molecules:
            log_file.write(molecule.log_entry())
    """

    # TODO: Revert original_ids to self.original_ids if we change the job
    #       monitoring to follow more object-oriented design. Otherwise,
    #       convert this to a static method or move the code to the
    #       execute method.
    def print_summary(self, original_ids, sample, elapsed_time=None):
        """Output a summary of the sample (number of molecules, time taken, etc.).

        :param sample: Molecule IDs from input molecule packet
        :param sample:
        :param elapsed_time:
        """
        if elapsed_time is not None:
            print(f"Step took {elapsed_time:.3} seconds")

        print(f"Sample has {len(sample)} molecules")
        parent_ids = set(str(m.molecule_id).split(".")[0] for m in sample)
        #TODO: this does not seem to be calculated correctly: always gives 0
        percent_original_represented = len(original_ids.intersection(parent_ids))/len(original_ids)
        print(f"Percent of the original ids that are still represented: {percent_original_represented:0.2%}")

        size_bin_cutoffs = [0, 100, 500, 1000, 1_000_000_000]
        size_counts, _ = np.histogram([len(mol) for mol in sample], bins=size_bin_cutoffs)
        print(f"Counts of molecules in size ranges:")
        for i in range(len(size_counts)-1):
            print(f" <={size_bin_cutoffs[i+1]}: {size_counts[i]}")
        print(f">{size_bin_cutoffs[-2]}: {size_counts[-1]}")


    def get_validation_attributes(self, configuration, global_config, directory_structure, input_molecule_packet):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the LibraryPrepPipeline job.

        NOTE: This method is included to satisfy the require requirements of the
              AbstractBeersPipeline interface, which this class implements. The output from
              this method is currently incompatible with the is_output_valid() method, and
              the validation attributes are assembled by the dispatcher, without any input
              from the class itself. This code is an intermediate step along the way to
              updating BEERS2 to fully support the BEERS_UTILS job submission framework. At
              the moment, the dispatcher uses a generic command to submit the library prep
              and sequencing pipelines, without instantiating objects for either pipeline.
              This means there's not currently any way to call this method from the
              dispatcher, since it's not a static method.
        NOTE: The is_output_valid() method currently requires the output_directory_path,
              which it uses to construct paths to the log and data directories. The log
              and data directories are currently stored as instance variables and would
              require a regex/strip to re-derive the output_directory_path.

        :param configuration: dictionary of the configuration data relevant to this pipeline stage.
                              [Note: this parameter is captured just so get_validation_attributes()
                              accepts the same arguments as get_commandline_call(). It is not used
                              here.]
        :param global_config: dictionary of full config data [Note: this parameter is captured just
                              so get_validation_attributes() accepts the same arguments as
                              get_commandline_call(). It is not used here.]
        :param directory_structure: instructions for navigating the subdirectories under the output
                                    directory
        :param input_molecule_packet: the molecule packet to run through this pipeline stage. Only
                                      the molecule ID is used.
        :return: Dictionary with a LibraryPrepPipeline's output_directory_path, molecule_packet_id,
                 and directory_structure.

        """

        validation_attributes = {}
        validation_attributes['data_directory_path'] = self.data_directory_path
        validation_attributes['log_directory_path'] = self.log_directory_path
        validation_attributes['packet_id'] = input_molecule_packet.molecule_packet_id
        validation_attributes['directory_structure'] = directory_structure
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of LibraryPrepPipeline for a specific job/execution is
        correctly formed and valid, given a job's output directory, molecule,
        packet ID, and directpry structure. Prepare these attributes for a given
        job using the get_validation_attributes() method.

        NOTE: This method is included to satisfy the require requirements of the
              AbstractBeersPipeline interface, which this class implements. This method is
              is currently incompatible with output from the get_validation_attributes()
              method. The validation attributes are assembled by the dispatcher, without any
              input from the class itself. This code is an intermediate step along the way
              to updating BEERS2 to fully support the BEERS_UTILS job submission framework.
              At the moment, the dispatcher uses a generic command to submit the library
              prep and sequencing pipelines, without instantiating objects for either
              pipeline. This means there's not currently any way to call the
              get_validation_attributes() method from the dispatcher, since it's not a
              static method.
        NOTE: The get_validation_attributes() method currently provides separate paths
              for the log and data directories, since it stores both of those paths as
              instance variables. This method derives both of these paths given the
              output_directory_path.

        Parameters
        ----------
        validation_attributes : dict
            A job's output_directory_path, molecule_packet_id, and directory_structure
            when running the library prep pipeline.

        Returns
        -------
        boolean
            True  - LibraryPrepPipeline output files were created and are well formed.
            False - LibraryPrepPipeline output files do not exist or are missing data.

        """

        output_directory_path = validation_attributes['output_directory_path']
        molecule_packet_id = validation_attributes['packet_id']
        directory_structure = validation_attributes['directory_structure']

        valid_output = False

        # Construct output filenames/paths
        log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        subdirectory_list = \
            GeneralUtils.get_output_subdirectories(molecule_packet_id, directory_structure)
        data_subdirectory_path = os.path.join(data_directory_path, *subdirectory_list)
        log_file_path = os.path.join(log_directory_path,
                                     LibraryPrepPipeline.pipeline_log_subdirectory_name,
                                     *subdirectory_list,
                                     f"{LibraryPrepPipeline.stage_name}_"
                                     f"molecule_pkt{molecule_packet_id}.log")
        results_filename = f"{LibraryPrepPipeline.stage_name}_" \
                           f"result_molecule_pkt{molecule_packet_id}.txt"
        results_file_path = os.path.join(data_subdirectory_path, results_filename)
        results_quant_filename = f"{LibraryPrepPipeline.stage_name}_" \
                                 f"result_molecule_pkt{molecule_packet_id}.quant_file"
        results_quant_file_path = os.path.join(data_subdirectory_path, results_quant_filename)

        if os.path.isfile(results_file_path) and \
           os.path.isfile(results_quant_file_path) and \
           os.path.isfile(log_file_path):

            #Read last line in log file
            line = ""
            with open(log_file_path, "r") as log_file:
                for line in log_file:
                    line = line.rstrip()
            if line == "Library prep pipeline completed successfully":
                valid_output = True

        return valid_output

    @staticmethod
    def main(seed, configuration, configuration_file_path, output_directory_path,
             directory_structure, molecule_packet_filename, packet_id,
             distribution_directory=None, molecules_per_packet_from_distribution=10000,
             sample_id=None,
             ):
        """
        This method would be called by a command line script in the bin directory.  It sets a random seed, loads a
        directory containing the relevant parts of the user's configuration file, unmarshalls a molecule packet from
        the provided molecule packet filename, initializes and validates the library prep pipeline stage and then
        executes it for the molecule packet.
        :param seed: value to use as the seed for the random number generator
        :param configuration: the json string containing the configration data specific to the library prep pipeline
        :param configuration_file_path: path to the full config file
        :param output_directory_path: top level output directory path for this pipeline stage
        :param directory_structure: instructions for creating the scaffolding needed to house the pipeline data and logs
        :param molecule_packet_filename: the file from which to unmarshall the molecule packet
        :param packet_id: id number to assign the packet
        :param distribution_directory: directory of CAMPAREE output distribution datas from which to generate molecules
            (default None, must supply molecule_packet_filename instead)
        :param molecules_per_packet_from_distribution: packet size to generate if using distributions
        :param sample_id: id number of the sample
        """

        print(f"Initializing with seed {seed}")
        np.random.seed(int(seed))
        configuration = json.loads(configuration)
        with open(configuration_file_path) as config_file:
            global_configuration = json.load(config_file)
        molecule_maker_parameters = global_configuration['molecule_maker_parameters']
        packet_id = None if packet_id == 'None' else int(packet_id)
        if molecule_packet_filename:
            molecule_packet = MoleculePacket.from_CAMPAREE_molecule_file(molecule_packet_filename, packet_id)
        elif distribution_directory:
            mm_log_path = output_directory_path # TODO should this go somewhere else?
            distribution_dir = pathlib.Path(distribution_directory)
            sample = Sample(
                sample_id = sample_id,
                sample_name = sample_id,
                adapter_sequences = [],
                fastq_file_paths = [],
                pooled = False,
            )
            molecule_maker = MoleculeMakerStep(
                    log_directory_path = mm_log_path,
                    parameters = molecule_maker_parameters,
            )
            molecule_packets = list(molecule_maker.execute(
                    sample = sample,
                    sample_data_directory = distribution_dir,
                    output_type =  "generator",
                    output_molecule_count = molecules_per_packet_from_distribution, # we only generate one packet
                    seed = seed,
                    molecules_per_packet = molecules_per_packet_from_distribution,
            ))
            molecule_packet = molecule_packets[0]
            molecule_packet.molecule_packet_id = packet_id
        else:
            raise BeersLibraryPrepValidationException("Neither molecule packet filename nor distribution directory provided")
        library_prep_pipeline = LibraryPrepPipeline(output_directory_path)
        library_prep_pipeline.execute(configuration, global_configuration, directory_structure, molecule_packet)

class BeersLibraryPrepValidationException(Exception):
    pass

if __name__ == '__main__':
    """
    Prepare and process command line arguments. This code is an intermediate
    step along to way to updating BEERS2 to full support the BEERS_UTILS job
    submission framework. Ideally, we can move this argument processing to
    the main function itself, and distribute some of the main() function's
    operations to the initialization and execute functions.
    """

    parser = argparse.ArgumentParser(description='Library Prep Pipeline')
    parser.add_argument('-s', '--seed', required=True, help="Seed", type=int)
    parser.add_argument('-c', '--config', required=True, help='Configuration')
    parser.add_argument('-C', '--config_file', required=True, help='Configuration File Path')
    parser.add_argument('-o', '--output_directory', required=True, help='Path to output directory.')
    parser.add_argument('-p', '--molecule_packet_filename', required=False, default=None, help="Serialized molecule packet filename (or none, if using distribution).")
    parser.add_argument('-d', '--directory_structure', required=True, help="Structure of data and logs directories.")
    parser.add_argument('--packet_id', required=True, help="ID of the packet to be used (either int or None if using id from the packet")
    parser.add_argument('-D', '--distribution_directory', required=False, default=None, help="Directory of CAMPAREE distributions (or none, if using molecule packet).")
    parser.add_argument('-N', '--num_molecules', required=False, default=10_000, help="Number of molecules to generate if using distributions.", type=int)
    parser.add_argument('-S', '--sample_id', required=False, default=None, help="Id number of the sample (used only when generating from distributions)", type=int)

    args = parser.parse_args()
    LibraryPrepPipeline.main(args.seed,
                             args.config,
                             args.config_file,
                             args.output_directory,
                             args.directory_structure,
                             args.molecule_packet_filename,
                             args.packet_id,
                             distribution_directory=args.distribution_directory,
                             molecules_per_packet_from_distribution=args.num_molecules,
                             sample_id=args.sample_id)
