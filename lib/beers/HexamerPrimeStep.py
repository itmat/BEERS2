from molecule import Molecule

class HexamerPrimeStep:

    def __init__(self, log_file, parameters):
        self.history_filename = log_file
        print("Hexamer_primer_step instantiated")

    def execute(self, sample):
        """
        Starting with a basic implementation that assumes all fragments get
        primed perfectly at their 3' ends.
        """
        print("Hexamer prime step starting")
        primed_sample = []
        primed_sample = sample
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in primed_sample:
                log_file.write(molecule.log_entry())
        print("Hexamer prime step complete")
        return primed_sample

    def validate(self):
        print("Hexamer prime step validating parameters")
        return True
