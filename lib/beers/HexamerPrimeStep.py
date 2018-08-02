class HexamerPrimeStep:

    def __init__(self):
        print("Hexamer_primer_step instantiated")

    def execute(self, sample):
        print(f"Hexamer prime step acting on sample {sample}")
        return f"Post Hexamer prime step sample"

    def validate(self, **kwargs):
        print(f"Hemamer prime step validating input kw args")
        return True