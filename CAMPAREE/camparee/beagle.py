import subprocess
import argparse
import sys
import os

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_utils import CampareeException
from camparee.camparee_constants import CAMPAREE_CONSTANTS

class BeagleStep(AbstractCampareeStep):

    #TODO: This is actually just the prefix of the output file. Beagle will generate
    #      two files in the data directory using this prefix: "beagle.vcf.log" and
    #      "beagle.vcf.vcf.gz". We should probably update the prefix so it at least
    #      doesn't contain the ".vcf" extension. However, before we do this, we need
    #      to make sure all downstream steps would use the updated filename.
    #TODO: Also, beagle creates its own log file. We should probably either suppress
    #      this or use beagle's own log file in place of the logging we currently do
    #      (which offers completely redundant information). Or, we update our own
    #      logging to track different things than beagle's native logging.
    BEAGLE_OUTPUT_FILENAME = CAMPAREE_CONSTANTS.BEAGLE_OUTPUT_PREFIX

    #Beagle takes input from the VariantsCompilationStep. This will make sure the
    #input filename matches up with the name used by the VariantsCompilationStep.
    BEAGLE_INTPUT_FILENAME = CAMPAREE_CONSTANTS.VARIANTS_COMPILATION_OUTPUT_FILENAME

    #Name of file where script logging stored.
    BEAGLE_LOG_FILENAME = CAMPAREE_CONSTANTS.BEAGLE_LOG_FILENAME

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path

    def validate(self):
        return True

    def execute(self, beagle_jar_path, seed=None):
        """
        Entry point into the beagle step. This ends up running the Beagle jar
        from the command line.

        Parameters
        ----------
        beagle_jar_path : string
            Path to the beagle JAR file.
        seed : int
            Seed for random number generator. Used so repeated runs will produce
            the same results.

        """
        input_file_path = os.path.join(self.data_directory_path, BeagleStep.BEAGLE_INTPUT_FILENAME)
        output_file_path = os.path.join(self.data_directory_path, BeagleStep.BEAGLE_OUTPUT_FILENAME)
        log_file_path = os.path.join(self.log_directory_path, BeagleStep.BEAGLE_LOG_FILENAME)
        command = f"java -jar {beagle_jar_path} gt={input_file_path} out={output_file_path}"
        if seed is not None:
            command += f" seed={seed}"
        with open(log_file_path, "w") as log_file:
            #TODO update the output here to make proper use of Python's logging
            #     module and functionality (it manages dual writing to both
            #     console and logging files).
            print(f"Calling beagle with command: {command}")
            print("Beagle output follows:\n")
            log_file.write(f"Calling beagle with command: {command}.\n\n")
            log_file.write("Beagle output follows:\n")

            try:
                beagle_result = subprocess.run(command, shell=True, check=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.STDOUT, #Redirect stderr to stdout.
                                               encoding="ascii")
            except subprocess.CalledProcessError as beagle_run_exception:
                log_file.write("\n*****ERROR: Beagle command failed:\n")
                log_file.write(f"\tExit code: {beagle_run_exception.returncode}\n")
                log_file.write("\n*****STDOUT:\n")
                log_file.write(f"{beagle_run_exception.stdout}\n")
                log_file.write("\n*****STDERR:\n")
                log_file.write(f"{beagle_run_exception.stderr}\n")
                raise CampareeException(f"\nBeagle process failed. "
                                        f"For full details see {log_file_path}\n")

            print(beagle_result.stdout)
            print(f"\nFinished running Beagle.\n")
            log_file.write(f"{beagle_result.stdout}\n")
            log_file.write(f"\nFinished running Beagle.\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, beagle_jar_path, seed=None):
        """
        Prepare command to execute the BeagleStep from the command line, given
        all of the arugments used to run the execute() function.

        Parameters
        ----------
        beagle_jar_path : string
            Path to the beagle JAR file.
        seed : int
            Seed for random number generator. Used so repeated runs will produce
            the same results.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the beagle.py script.
        beagle_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        beagle_step_path = beagle_step_path.rstrip('c')

        command = (f" python {beagle_step_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --beagle_jar_path {beagle_jar_path}")

        if seed is not None:
            command += f" --seed {seed}"

        return command

    def get_validation_attributes(self, beagle_jar_path, seed=None):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the BeagleStep job.

        Parameters
        ----------
        beagle_jar_path : string
            Path to the beagle JAR file. [Note: this parameter is captured just
            so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]
        seed : int
            Seed for random number generator. Used so repeated runs will produce
            the same results. [Note: this parameter is captured just so
            get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A BeagleStep run's data_directory and log_directory.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        return validation_attributes

    @staticmethod
    def main():
        """
        Entry point into script. Allows script to be executed/submitted via the
        command line.
        """

        parser = argparse.ArgumentParser(description='Command line wrapper around'
                                                     ' the Beagle step')
        parser.add_argument('--log_directory_path')
        parser.add_argument('--data_directory_path')
        parser.add_argument('--beagle_jar_path')
        parser.add_argument('--seed', type=int, default=None)
        args = parser.parse_args()

        beagle_step = BeagleStep(args.log_directory_path,
                                 args.data_directory_path)
        beagle_step.execute(beagle_jar_path=args.beagle_jar_path,
                            seed=args.seed)

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of BeagleStep for a specific job/execution is correctly
        formed and valid, given the run's data and log directories. Prepare these
        attributes using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A CAMPAREE run's data_directory and log_directory.

        Returns
        -------
        boolean
            True  - BeagleStep output files were created and are well formed.
            False - BeagleStep output files do not exist or are missing data.

        """
        data_directory = validation_attributes['data_directory']
        log_directory = validation_attributes['log_directory']

        valid_output = False

        #The way beagle is set to run above, it generates an output file based
        #on the given filename prefix, with a ".vcf.gz" suffix added.
        output_file_path = os.path.join(data_directory, BeagleStep.BEAGLE_OUTPUT_FILENAME + ".vcf.gz")
        log_file_path = os.path.join(log_directory, BeagleStep.BEAGLE_LOG_FILENAME)
        if os.path.isfile(output_file_path) and os.path.isfile(log_file_path):
            #Read last line in beagle log file
            line = ""
            with open(log_file_path, "r") as log_file:
                for line in log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                valid_output = True

        return valid_output

if __name__ == "__main__":
    sys.exit(BeagleStep.main())
