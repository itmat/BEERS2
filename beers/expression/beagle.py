import subprocess
import argparse
from timeit import default_timer as timer
import sys
import os

class BeagleStep:

    def __init__(self, logfile, data_directory_path, parameters):
        self.data_directory_path = data_directory_path

    def validate(self):
        return True

    def execute(self, resources):
        third_party_software_directory = os.path.join(resources['resources_folder'], "third_party_software")
        beagle_filename = [filename for filename in os.listdir(third_party_software_directory)
                           if "beagle" in filename][0]
        beagle_file_path = os.path.join(third_party_software_directory, beagle_filename)
        input_file_path = os.path.join(self.data_directory_path, "all_variants.vcf")
        output_file_path = os.path.join(self.data_directory_path, "beagle.vcf")
        command = f"java -jar {beagle_file_path} gt={input_file_path} out={output_file_path}"
        print(f"Calling beagle with command: {command}")
        result = subprocess.call(command, shell=True)
        print(f"Finished running Beagle.")
        return result
