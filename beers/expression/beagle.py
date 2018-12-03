import subprocess
import argparse
from timeit import default_timer as timer
import sys

class BeagleStep:

    def __init__(self, logfile, data_directory_path, parameters):
        self.command = parameters["command"]

    def validate(self):
        return True

    def execute(self):
        print("Calling BEAGLE")
        result = subprocess.call(
            f"java -classpath {self.command}", shell=True)
        return result

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Beagle')
        required_named = parser.add_argument_group('required named arguments')
        required_named.add_argument('-c', '--command', required=True, help='Command line needed to run Beagle')
        args = parser.parse_args()
        print(args)

        parameters = {
            'command': args.command
        }

        beagle_step = BeagleStep(parameters)
        start = timer()
        beagle_step.execute()
        end = timer()
        sys.stderr.write(f"Beagle Step: {end - start} sec\n")

if __name__ == "__main__":
    sys.exit(BeagleStep.main())


'''Example call
python beagle.py \
 -c "java -classpath /home/crislawrence/eclipse-workspace/Stubs/bin org.crislawrence.stubs.MockProcess 5"
'''