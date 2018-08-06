class HexamerPrimeStep:

    def __init__(self, log_file, parameters):
        print("Hexamer_primer_step instantiated")

    def execute(self, sample):
        print("Hexamer prime step acting on sample")
        return f"Post Hexamer prime step sample"

    def validate(self):
        print(f"Hexamer prime step validating parameters")
        return True
