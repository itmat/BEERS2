import subprocess
import argparse
from timeit import default_timer as timer
import sys
import os

class BeagleStep:

    def __init__(self, logfile, data_directory_path, parameters=dict()):
        self.data_directory_path = data_directory_path

    def validate(self):
        return True

    def execute(self, beagle_file_path):
        input_file_path = os.path.join(self.data_directory_path, "all_variants.vcf")
        output_file_path = os.path.join(self.data_directory_path, "beagle.vcf")
        command = f"java -jar {beagle_file_path} gt={input_file_path} out={output_file_path}"
        print(f"Calling beagle with command: {command}")
        result = subprocess.call(command, shell=True)
        print(f"Finished running Beagle.")
        return result
